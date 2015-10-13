% GSSM  General state space model for camera gps imu integration
%
%   The following state space model is used :
% the state vector is defined as delta(ecef x,y,z), delta(vn, ve, vd),
% phi (roll, pitch, yaw), accelerometer bias and gyro bias), for srcdkf
%   and ekf. 
% attitude expressed in quaternion, since accelerometer and gyro scale errors are
% relatively small, only bias are accounrted for and modeled as a
%   first order Gaussian Markov since the driving noise covariance is
%   unattainable given the bias specs of the sensors
% note each particle state has 16 dimensions because of quaternions
%   And the observations at time k, O(1,2)(k) are :
%           1, the GPS derived ecef position
%           2,- the camera matched point pairs
%
%   The state dynamics are driven by a 12 (accel and gyro white noise and bias white noise)
%   dimensional white Gaussian noise source and the
%   observations are corrupted by additive white Gaussian noise.
%=============================================================================================

function [varargout] = model_interface(func, varargin)

  switch func

    %--- Initialize GSSM data structure --------------------------------------------------------
    case 'init'
      model = init(varargin{1});
        error(consistent(model,'gssm'));               % check consistentency of initialized model
      varargout{1} = model;

    %--------------------------------------------------------------------------------------------
    otherwise

      error(['Function ''' func ''' not supported.']);

  end


%===============================================================================================
function model = init(options)
 
  model.type = 'gssm';                                      % object type = generalized state space model
  model.tag  = 'GSSM_Cam_IMU_GPS_Tracking';  % ID tag
  model.InvalidateIMUerrors=options.InvalidateIMUerrors;
  model.imutype=options.imutype;
  if (options.useCam)
  model.camtype=options.camtype;
  model.camLA=options.camLA;%lever arm offset, the coordinates of the lens in the IMU frame
  model.Rb2c=[1 0 0; 0 0 1; 0 -1 0]*[0 1 0; -1 0 0; 0, 0, 1];
  end
  model.LA=options.LA;% lever arm offset, the coordinates of the GPS antenna in the IMU frame
   
  model.statedim   = 16;                      %   state dimension, 15+1 for quarternion
  model.obsdim     = 3;                      %   observation dimension, not yet used in pf and sppf
  model.U1dim      = 8;                      %   exogenous control input 6 dimension,
  % imu accel and angular rate, NOT delta v and delta theta, previous
  % epoch, current epoch time
  model.U2dim      = 1;                      %   exogenous control input 1 dimension, detemine which hfun to use
  % U2=-1 NO observation, 1 GPS, 2 ZUPT, 3, NHC, 4 Camera
  model.Vdim       = 12;                     %   process noise dimension, 
  %accel and gyro white noise and accel and gyro bias white noise
  model.oModelNr=4; %number of observation models
 
  model.Ndim       = [3;3;2;1];             %   observation noise dimension for GPS, ZUPT, NHC and camera 
  model.h=sqrt(3);    % scale factor (CDKF parameter h)
  
  model.ffun      = @ffun;                   % file handle to FFUN
  model.hfun      = @hfun;                   % file handle to HFUN
  model.prior     = @prior;
  model.likelihood = @likelihood;            % file handle to LIKELIHOOD
  model.innovation = @innovation;            % file handle to INNOVATION
  model.generateparticles=@generateparticles;
  model.estimation=@estimation;
  model.diffstate=@diffstate;
  model.hfuncam=@hfuncam;
  model.measurementupdate=@measurementupdate;
  model.MU_AdditiveNoise=@MU_AdditiveNoise;
  model.setparams  = @setparams;             % file handle to SETPARAMS
  model.paramdim   = 0;   
  model.params = 0;
  
  Arg.type = 'gaussian';
  Arg.dim = model.Vdim;
  Arg.mu  = zeros(Arg.dim,1);
  Arg.cov_type = 'full';
  IMU_ERRDEF=imu_err_defs_v000(options.imutype);
  model.IMU_ERRDEF=IMU_ERRDEF;
%   beta=1e-4;
  Arg.cov= diag([IMU_ERRDEF.acc_vrw, IMU_ERRDEF.acc_vrw, IMU_ERRDEF.acc_vrw, ... % for accel white noise
      IMU_ERRDEF.gyro_arw, IMU_ERRDEF.gyro_arw, IMU_ERRDEF.gyro_arw, ... % for gyro white noise    
      [IMU_ERRDEF.acc_bias_Q, IMU_ERRDEF.acc_bias_Q, IMU_ERRDEF.acc_bias_Q],...% for accel random walk bias noise 
 [IMU_ERRDEF.gyro_bias_Q, IMU_ERRDEF.gyro_bias_Q, IMU_ERRDEF.gyro_bias_Q]])*options.dt;% for gyro random walk bias noise
%  1e-4 kind of like 2*beta with beta=1/correlation time constant
% 1e-3 and 1e-2 gives almost the same results for Personal navigator MEMS,
% 3DM GX 3-35, but 1e-4 seems the best, and 1 seems the worst.
 model.pNoise = gennoiseds(Arg);            % process noise : zero mean white Gaussian noise, cov = (1e-3)^2(dynamics) and (1e-4)^2 (tone frequency, very stable)

P0 =zeros(model.statedim-1);  % Covariance for the initial state
P0(1:3,1:3)=eye(3)*5^2; %position errors in meter, position in e frame
P0(4:6,4:6)=eye(3)*0.2^2;
P0(7:9,7:9)=diag([options.initAttVar,options.initAttVar,options.initAttVar*2].^2); 
%attitude errors in radian (may be determined later during coarse alignment)
P0(10:12,10:12)=eye(3)*IMU_ERRDEF.initacc_bias_err^2;
P0(13:15,13:15)=eye(3)*IMU_ERRDEF.initgyro_bias_err^2;
model.Sx=chol(P0)';

  %this noise for sigma point particle filter to compute prior
  Arg.type = 'gaussian';
  Arg.dim = model.statedim-1; % -1 because of the quaternion representing attitude
  Arg.mu  = zeros(Arg.dim,1);
  Arg.cov_type = 'full';
  Arg.cov= eye(Arg.dim);
  model.pNoiseSP = gennoiseds(Arg);   
  
model.oNoise = cell(model.oModelNr,1); 
for i=1:size(model.Ndim,1)
  Arg.type = 'gaussian';
  Arg.dim = model.Ndim(i);
  Arg.mu = zeros(Arg.dim,1);
  Arg.cov_type ='full';
  switch i
      case 1%gps
          Arg.cov  =diag([0.05, 0.05,0.1].^2);  % to be changed for each GPS measurement
      case 2% zupt
          Arg.cov  =eye(3)*options.sigmaZUPT^2;
      case 3% NHC
          Arg.cov  =eye(2)*options.sigmaNHC^2;
      case 4% camera
          Arg.cov  =eye(1)*options.sigmaCAM^2;
      otherwise
          error([' Unknown measurement noise type in gssm_cam_ins_gps.m']);
  end
  model.oNoise{i} = gennoiseds(Arg);            % observation noise : zero mean white Gaussian noise, cov=0.0175^2 (bearings=1° error) and 0.06^2 (frequencies)
end
function model = setparams(model, params, index_vector)
model.params=params;
% dummy function
%===============================================================================================
function particles = generateparticles(model,x0, NrPart)
% state vector contains ecef xyz, vn, ve,vd, roll pitch yaw, accelerometer and gyro bias
switch model.inftype
    case {'sppf','pf'}
        %to be corrected
%         particles=zeros(model.statedim,NrPart);%initialize particles
%         % (1) Initialize
%         xnp= repmat(x0,1,NrPart) + model.Sx*randn(model.statedim-1,NrPart);  % Nonlinear states
%         
%         for counter=1:NrPart
%             particles(:,counter)=[xnp(1:6,counter);att2qua(xnp(7:9,counter));xnp(10:end,counter)];
%         end
    case {'srcdkf','srukf'}
        %remember to keep model.Sx updated before calling this function
        nsp1 = 2*(model.statedim-1+model.Vdim)+ 1;          % number of sigma points
        Z   = cvecrep([x0; model.pNoise.mu],nsp1);
        Xdim=model.statedim-1;
        Vdim=model.Vdim;
        Zeros_Xdim_X_Vdim = zeros(Xdim,Vdim);
        Zeros_Vdim_X_Xdim = zeros(Vdim,Xdim);
        Sv =  chol(model.pNoise.cov)';         % matrix square root of process noise covariance
        Sz  = [model.Sx Zeros_Xdim_X_Vdim; Zeros_Vdim_X_Xdim Sv];
        hSz = model.h*Sz;
        hSzM = [hSz -hSz];
        for counter=2:nsp1
            Z(7:10,counter)=quatmult_v001(rvec2quat_v000(hSzM(7:9,counter-1)), Z(7:10,counter),0);
            Z(1:6,counter)=Z(1:6,counter)+hSzM(1:6,counter-1);
            Z(11:end,counter)=Z(11:end,counter)+hSzM(10:end,counter-1);            
        end
        particles=Z;
end

%===============================================================================================
%-- State transition function (vehicle dynamic model)
% with accelerometer and gyro input at time k-1 to predict state at k, 
% i.e., x(k|k-1) from state at k-1, U1 contains accelerometer measurement
% acceleration, not delta v, and gyro output angular rate, not delta Theta
% and 7th row is the previous epoch and 8th is the current time
% V is composed of wa, wg, wab, wgb in each column
function new_state = ffun(model, state, V, U1)
NrPart=size(state,2);
new_state=zeros(size(state));
dt=U1(end)-U1(end-1);
if isempty(V) % deterministic propogation
    for i = 1:NrPart
        gyroinc=U1(4:6,i)*dt;
        accinc=U1(1:3,i)*dt;
        angleinc=gyroinc;
        vr_a=0.5*cross(gyroinc, accinc);
        velinc=accinc+vr_a;
        %Calibrate the imu
        if(model.InvalidateIMUerrors)
            spuracc=find(abs(state(11:13,i))>model.IMU_ERRDEF.initacc_bias_err,1);
            spurgyro=find(abs(state(14:16,i))>model.IMU_ERRDEF.initgyro_bias_err,1);
            if(~isempty(spuracc)||~isempty(spurgyro))
                disp(['IMU errors diverge at ' num2str(U1(end)) ' particle ' num2str(i) '!']);
                state(11:13,i)=0;
                state(14:16,i)=0;
            end
        end        
        velinc=velinc-state(11:13,i)*dt;
        angleinc=angleinc-state(14:16,i)*dt;
        
        %run strapdown with angle and velocity increments
        [qbn_new, Vn_new, llh_new]=strapdown_ned_quat_v000(state(7:10,i), state(4:6,i),ecef2geo_v000(state(1:3,i),0),velinc, angleinc, dt);
        new_state(1:10,i)=[ecef2geo_v000(llh_new,1);Vn_new;qbn_new];      
        new_state(11:end,i)=state(11:end,i);
    end
else
    for i = 1:NrPart      
        gyroinc=U1(4:6,i)*dt;
        accinc=U1(1:3,i)*dt;
        angleinc=gyroinc;
        vr_a=0.5*cross(gyroinc, accinc);
        velinc=accinc+vr_a;
        %Calibrate the imu
        if(model.InvalidateIMUerrors)
            spuracc=find(abs(state(11:13,i))>model.IMU_ERRDEF.initacc_bias_err,1);
            spurgyro=find(abs(state(14:16,i))>model.IMU_ERRDEF.initgyro_bias_err,1);
            if(~isempty(spuracc)||~isempty(spurgyro))
                disp(['IMU errors diverge at ' num2str(U1(end)) ' particle ' num2str(i) '!']);
                state(11:13,i)=0;
                state(14:16,i)=0;
            end
        end
        
        velinc=velinc-(state(11:13,i)*dt+V(1:3,i));
        angleinc=angleinc-(state(14:16,i)*dt+V(4:6,i));
        
        %run strapdown with angle and velocity increments
       [qbn_new, Vn_new, llh_new]=strapdown_ned_quat_v000(state(7:10,i), state(4:6,i),ecef2geo_v000(state(1:3,i),0),velinc, angleinc, dt);
        new_state(1:10,i)=[ecef2geo_v000(llh_new,1);Vn_new;qbn_new];     
        new_state(11:end,i)=state(11:end,i)+V(7:end,i)-[state(11:13,i)/model.IMU_ERRDEF.acc_bias_Tc;state(14:end,i)/model.IMU_ERRDEF.gyro_bias_Tc]*dt;
    end
end
%compute the difference between two states in terms of ecef coord, delta
%Vned, delta attitude and delta the following fields X1-X2
% X1 and X2 are composed of ecef xyz, vned, quaternion and so on
function delta=diffstate(model, X1, X2)
NrPart=size(X1,2);
delta=zeros(model.statedim-1,NrPart);
for i=1:NrPart
    delta(7:9,i)=quat2rot_v000(quatmult_v001(X1(7:10,i),X2(7:10,i),2));% kind of like subangle(qua2att(X1(7:10,i)),qua2att(X2(7:10,i)))
end
delta(1:6,:)=X1(1:6,:)-X2(1:6,:);
delta(10:end,:)=X1(11:end,:)-X2(11:end,:);

%===============================================================================================
%-- State observation function
% state of n particles in each column, predict the observations k from state
% at k
function observ = hfun(model, state, N, U2)
NrPart=size(state,2);
switch U2(1)
    case 1
        observ_ = zeros(3,NrPart); % in ecef frame
        for i=1:NrPart
            imugeo=ecef2geo_v000(state(1:3,i),0);
            Cen=llh2dcm_v000(imugeo(1:2,1),[0,1]);
            observ_(:,i)=state(1:3,i)+Cen'*quatrot_v000(state(7:10,i),model.LA,0);
        end
    case 2% zupt
        observ_ = state(4:6,:); % in ecef frame
    case 3%nhc
        observ_ = zeros(2,NrPart); % in body frame, right and down velocity
        for i=1:NrPart
            %                        Cnb=quat2dcm_v000(particlesPred(7:10,i))';
            %                        OBS(:,i)=Cnb(2:3,:)*particlesPred(4:6,i);
            fool=quatrot_v000(state(7:10,i),state(4:6,i),1);
            observ_(:,i)=fool(2:3);
        end        
    otherwise   
        error(' Unknown measurement type in hfun() gssm_cam_ins_gps.m');
  
end
if isempty(N),
    observ = observ_;
else
      observ = observ_ + N;
end




function observ = hfuncam(model, state, lastState, camdata)
NrPart=size(state,2);
observ=zeros(size(camdata,2),NrPart);
%predict the observation, current frame, C2, previous frame, C1,
imugeo=ecef2geo_v000(lastState(2:4),0);
Cen=llh2dcm_v000(imugeo(1:2,1),[0,1]);
xyz_lenspre=lastState(2:4)+Cen'*quatrot_v000(lastState(8:11),model.camLA,0);
Cnb=quat2dcm_v000(lastState(8:11))';
re2c1=model.Rb2c*Cnb*Cen;
cam = initialize_cam(model.camtype);
for num=1:NrPart
    imugeo=ecef2geo_v000(state(1:3, num),0);
    Cen=llh2dcm_v000(imugeo(1:2,1),[0,1]);
    xyz_lenscur=state(1:3, num)+Cen'*quatrot_v000(state(7:10, num),model.camLA,0);
    Cnb=quat2dcm_v000(state(7:10, num))';
    re2c2=model.Rb2c*Cnb*Cen;
%     tcur2pre=re2c1*(xyz_lenscur-xyz_lenspre);
    tcur2pre=model.Rb2c*quatrot_v000(lastState(8:11),Cen*(xyz_lenscur-xyz_lenspre),1)
    tcur2pre=tcur2pre/norm(tcur2pre,2);
    rcur2pre=re2c1*re2c2';
    % compute thepredicted observations for each pair    
    for colt=1:size(camdata,2)
        nickel=cam.K\[camdata(6:7,colt);1];
%         nickel= undistor_a_point( nickel, cam );   
        camdata(6:7,colt)=nickel(1:2);       
        nickel=cam.K\[camdata(9:10,colt);1];
        camdata(9:10,colt)=nickel(1:2);     
        observ(colt,num)=[camdata(6:7, colt)',1]*skew(tcur2pre)*rcur2pre*[camdata(9:10, colt);1];
    end    
end
%===============================================================================================
function tranprior = prior(model, V, U1, pNoiseDS)
tranprior = pNoiseDS.likelihood( pNoiseDS, V);

%===============================================================================================
function likeld = likelihood(model, obs, state, U2, oNoiseDS)

observ =  hfun(model, state, [], U2);
N = obs - observ;
% Calculate log likelihood
likeld = oNoiseDS{U2(1)}.likelihood( oNoiseDS{U2(1)}, N);

%===============================================================================================
function innov = innovation(obs, observ)
innov = obs - observ;
%=============================================================================================
function estimate=estimation(model, particles, weights)
if(size(weights,1)==2)
    tempw=ones(size(particles,2),1)*weights(2);
    tempw(1)=weights(1);
    weights=tempw;
end
%OUTPUT ESTIMATED state vector, ecef xyz, vned, rpy, etc
estimate= sum(rvecrep(weights',model.statedim).*particles,2);          % expected mean
% estimate roll picch yaw
estimate(7:10,1)=weighted_avg_quat(particles(7:10,:)',weights);
%input: Sx_, the cholesky square root of the state covriance before the
%measurement received, in matlab Sx_=chol(P0)', OBS, observation column
%vector, W1, the weights set 2 for averaging sigma points, U2, measurement
%type, Xh, the last predicted state X(k|k-1)
function  [Xh_new,Sx_new]=measurementupdate(model,Xh, Sx_,OBS,W2,U2)
Xdim=model.statedim-1;
Ndim=model.Ndim(U2);
nsp2   = 2*(Xdim+Ndim) + 1;          % number of sigma points (second set)
Zeros_Xdim_X_Ndim = zeros(Xdim,Ndim);
Zeros_Ndim_X_Xdim = zeros(Ndim,Xdim);
Sn=  chol(model.oNoise{U2}.cov)';
Z  = cvecrep([Xh; model.oNoise{U2}.mu] ,nsp2);
Sz = [Sx_ Zeros_Xdim_X_Ndim; Zeros_Ndim_X_Xdim Sn];
hSz = model.h*Sz;
hSzM = [hSz -hSz];
for counter=2:nsp2
    
    Z(7:10,counter)=quatmult_v001(rvec2quat_v000(hSzM(7:9,counter-1)), Z(7:10,counter),0);
    Z(1:6,counter)=Z(1:6,counter)+hSzM(1:6,counter-1);
    Z(11:end,counter)=Z(11:end,counter)+hSzM(10:end,counter-1);
    
end
% build sigma-point set

%-- Calculate predicted observation mean
Y_ = model.hfun( model, Z(1:Xdim+1,:), Z(Xdim+2:Xdim+Ndim+1,:), U2);
yh_ = W2(1,1)*Y_(:,1) + W2(1,2)*sum(Y_(:,2:nsp2),2);
C = W2(2,1) * ( Y_(:,2:Xdim+Ndim+1) - Y_(:,Xdim+Ndim+2:nsp2) );
D = W2(2,2) * ( Y_(:,2:Xdim+Ndim+1) + Y_(:,Xdim+Ndim+2:nsp2) - cvecrep(2*Y_(:,1),Xdim+Ndim));
[temp,Sy] = qr([C D]',0);
Sy = Sy';
% MEASUREMENT UPDATE

Syx1 = C(:,1:Xdim);
Syw1 = C(:,Xdim+1:end);

Pxy = Sx_*Syx1';

KG = (Pxy / Sy') / Sy;
inov =OBS-yh_;  % this coincides with Shin, 2005
upd = KG*inov;
% becuase we estimate the error state of re, vn and
% rpy in srcdkf, we need to subtract the -upd from the estimated
% values of them, and we estimate bias and scale factors directly in
% srcdkf, so we just add the -upd to the state vector.
% THEREFORE, we need to distinguish the Xh(1:10) and
% Xh(11:end), Xh(11:end) is in the state vector, Xh(1:10)
% is not, but it turns out Shin, 2005 was correct,we need to subtract -upd
%from the bias vector. I still need to investigate more about this
Xh_new=Xh;
Xh_new(1:6)=Xh(1:6)+upd(1:6);
Xh_new(7:10)=quatmult_v001(rvec2quat_v000(upd(7:9)),Xh(7:10),0);
Xh_new(11:end)=Xh(11:end)+upd(10:end);

[temp,Sx] = qr([Sx_-KG*Syx1 KG*Syw1 KG*D]',0);
Sx_new=Sx';
%measurement update for additive noise, for the camera case, U2 is the cam
%point pairs, camdata contains the uncalibrated unnormalized camera points
%coordinates in U,V frame with upper left as the origin
function  [Xh_new,Sx_new]=MU_AdditiveNoise(model,Xh, Sx_,OBS,W2,lastState, camdata)
U2=4;
Xdim=model.statedim-1;
nsp2   = 2*Xdim+1;          % number of sigma points (second set)
numPairs=size(camdata,2);
Sn=diag(repmat(sqrt(model.oNoise{U2}.cov),numPairs,1));
% build sigma-point set
Z  = cvecrep(Xh ,nsp2);
hSz = model.h*Sx_;
hSzM = [hSz -hSz];
for counter=2:nsp2    
    Z(7:10,counter)=quatmult_v001(rvec2quat_v000(hSzM(7:9,counter-1)), Z(7:10,counter),0);
    Z(1:6,counter)=Z(1:6,counter)+hSzM(1:6,counter-1);
    Z(11:end,counter)=Z(11:end,counter)+hSzM(10:end,counter-1);    
end
%-- Calculate predicted observation mean
Y_ = model.hfuncam( model, Z(1:Xdim+1,:),lastState, camdata);
yh_ = W2(1,1)*Y_(:,1)+ W2(1,2)*sum(Y_(:,2:nsp2),2);
C = W2(2,1) * ( Y_(:,2:Xdim+1) - Y_(:,Xdim+2:nsp2) );
D = W2(2,2) * ( Y_(:,2:Xdim+1) + Y_(:,Xdim+2:nsp2) - cvecrep(2*Y_(:,1),Xdim));
[temp,Sy] = qr([C D Sn]',0);
Sy = Sy';
% MEASUREMENT UPDATE

Syx1 = C(:,1:Xdim);
Syw1 = C(:,Xdim+1:end);

Pxy = Sx_*Syx1';

KG = (Pxy / Sy') / Sy;
inov =OBS-yh_;  % this coincides with Shin, 2005
upd = KG*inov;
% becuase we estimate the error state of re, vn and
% rpy in srcdkf, we need to subtract the -upd from the estimated
% values of them, and we estimate bias and scale factors directly in
% srcdkf, so we just add the -upd to the state vector.
% THEREFORE, we need to distinguish the Xh(1:10) and
% Xh(11:end), Xh(11:end) is in the state vector, Xh(1:10)
% is not, but it turns out Shin, 2005 was correct,we need to subtract -upd
%from the bias vector. I still need to investigate more about this
Xh_new=Xh;
Xh_new(1:6)=Xh(1:6)+upd(1:6);
Xh_new(7:10)=quatmult_v001(rvec2quat_v000(upd(7:9)),Xh(7:10),0);
Xh_new(11:end)=Xh(11:end)+upd(10:end);
[temp,Sx] = qr([Sx_-KG*Syx1 KG*Syw1 KG*D]',0);
Sx_new=Sx';

                
