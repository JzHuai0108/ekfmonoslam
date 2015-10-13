% A better version of AddIMUErr_v001 (I will depreciate the _v001 version.)

% This script assumes that continuous Output of each sensor is modelled as:
%   x_dot=A*x+Bw
%   y=u+Cti*x+(M*u_all)'*Ctv*x+Dv
%   where,
%       x is the error states, e.g., $[b_a, b_g, s_a, s_g]$, $[b_a, b_g,
%       r_{a1}, r_{a1}, r_{a2}, r_{g1}, r_{g1}, r_{g2}]$ where s_g
%       represents the scale factor error of three gyros, r_{a2} represents
%       the misalignment of accelerometer in the y axis
%       E[|Bw|^2]=Q\delta(t-t'), E[|Bwdt|^2]=Qdt,
%       Cti extracts relevant bias error terms from x
%       y is a component of the sensor output, e.g., $[\tilde{a}_s^s, \tilde{\omega}_{is}^s]$
%       u is a component of true acceleration/rotation rate input
%       u_all=true_inp(:,in) all true input values at an epoch
%       $E[|Dv|^2]=R\delta(t-t')$, note Dv is continuous noise with noise
%       density, D^2=R, given on the specs of IMU. In generating IMU noise,
%       we have to discretize Dv, let Disv=\int_{0}^{\Delta
%       t}Dvdt/\Delta t, so E[|Disv|^2]=R/\Delta t
%       M extracts relevant variables from input vector in computing the
%       scale/misalignment effect, Ctv extracts relevant scale factor
%       error/misalignment error terms from x
%       The deafult value of M*u_all is u, where u_all is an input vector.

% true_inp, 6xN matrix, six rows represent $[a_s^s, \omega_{is}^s]$
% ErrDefs is an array of cell. Each cell defines the error model of the corresponding row of true_inp.
% Each member of ErrDefs must be as follows:
%   ErrDefs{i}.Cti, e.g., [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
%   ErrDefs{i}.Ctv  (This can be undefined or empty if there is no scale
%   factor or cross axis errors), e.g., [0,0,0, 0,0,0, 1, 0,0, 0,0,0]
%   ErrDefs{i}.R, where D*D'=R
%   ErrDefs{i}.M (M can be undefined or empty. Unless you use SRIMU,
%   completely ignore M.), e.g., [zeros(3), eye(3)]
%   The number of struct and number of rows of true_inp MUST be the same.
% StateErrDef, for all error terms in x,
%   StateErrDef.A
%   StateErrDef.Q, where B*B'=Q in continuous case
% STD0=sqrt(cov(x_0))
% dt is the sampling period of true_inp.

% sensor_out, 6xN matrix, its rows represent $[\tilde{a}_s^s, \tilde{\omega}_{is}^s]$
% sensor_err, kxN matrix, its rows corresponds to IMU errors

% This script ignores squaring type of errors (2nd order effects). Perhaps, I may add these effects
%  in a future version.
% This script also is not capable of simulating temperature dependent errors. (_v001 version was
%  capable of doing that.)

function [sensor_out, sensor_err]= AddIMUErr(true_inp, ErrDefs, StateModelDef, STD0, dt)
%Check the arguments' integrity
[nsen, ndat]=size(true_inp); %Each row corresponds a different sensor
nsta= size(STD0,1);   %total number of states
if (nsen~=length(ErrDefs))
    disp('The number of cells and number of rows of true_inp must be the same');
    return;
end
%Generate M matrices unless they are defined or not empty.
%The default M matrix corresponds to scale factor error
scale_err=zeros(nsen,1);
for in=1:nsen
    if isfield(ErrDefs{in}, 'Ctv') && ~isempty(ErrDefs{in}.Ctv)
        scale_err(in)=1;
        %Define the default value of M
        if (~isfield(ErrDefs{in}, 'M')) || isempty(ErrDefs{in}.M)
            ErrDefs{in}.M=zeros(1,size(true_inp,1));
            ErrDefs{in}.M(:,in)=1;
        end
    else
        scale_err(in)=0;
    end
end

%Convert to discrete time
%the model is in continuous time. Convert it to discrete time
[StateModelDef.A, StateModelDef.Q, ~]=dc2dc_v000(StateModelDef.A, StateModelDef.Q, [],dt,0,[]);
for in=1:nsen
    ErrDefs{in}.R=ErrDefs{in}.R/dt;
end
%Convert covariances into coefficients and define M matrices
[L,D]=ldl(StateModelDef.Q);
StateModelDef.B=L*D.^0.5;
for in=1:nsen
    %Coefficients for IMU error state generation
    [L,D]=ldl(ErrDefs{in}.R);
    ErrDefs{in}.D=L*D.^0.5;
end
%Output variables
sensor_out=zeros(nsen,ndat);
sensor_err=zeros(nsta,ndat);
%Generate and add errors
%Combine all models into a single matrix
A=StateModelDef.A;
B=StateModelDef.B;
Cti=zeros(nsen, nsta);
Ctv_scaled=zeros(nsen, nsta);
D=[];
for in=1:nsen
    Cti(in, :)=ErrDefs{in}.Cti;
    if (scale_err(in))
        Ctv_scaled(in,:)=(ErrDefs{in}.M*true_inp(:,1))'*(ErrDefs{in}.Ctv);
    else
        Ctv_scaled(in,:)=zeros(1,nsta);
    end
    D=diagmat_v000(ErrDefs{in}.D,D,1);
end

sensor_err(:,1)=STD0*randn(nsta,1);
sensor_out(:,1)=true_inp(:,1)+(Cti+Ctv_scaled)*sensor_err(:,1)+D*randn(nsen,1);
%Generate the rest
for in=2:ndat
    %Time varying component
    for k=1:nsen
        if (scale_err(k))
            Ctv_scaled(k, :)=(ErrDefs{k}.M*true_inp(:,in))'*(ErrDefs{k}.Ctv);
        else
            Ctv_scaled(k,:)=zeros(1,nsta);
        end
    end
    %Generate and add errors
    sensor_err(:,in)=A*sensor_err(:,in-1)+B*randn(nsta,1);
    sensor_out(:,in)=true_inp(:,in)+(Cti+Ctv_scaled)*sensor_err(:,in)+D*randn(nsen,1);
end