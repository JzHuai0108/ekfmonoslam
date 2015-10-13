%allan_std_par=[tao at peak, std value at peak;...];
%fs:discretization frequency used for A and B
function [A,B,P,allan_std]=allan2markov(allan_std_par, fs, t_axis)
nmar=size(allan_std_par,1); %number of markov to be computed
A=zeros(nmar,1);  
B=zeros(nmar,1);
P=zeros(nmar,1);
allan_std=zeros(nmar,length(t_axis));
for in=1:nmar
    %from allan par to markov par
    Tc=allan_std_par(in,1)/1.89;
    qc=allan_std_par(in,2)/0.437/sqrt(Tc);
    
    if (~isempty(fs))
        %convert the results into discrete form
        [cm dm]=markov1st_v000(-1/Tc,qc,[],1,fs);
        A(in)=dm(1);
        B(in)=dm(2);
        P(in)=dm(3);
    else
        %%Results in continuous form
        A(in)=-1/Tc;
        B(in)=qc;
        P(in)=qc^2/(-2*(-1/Tc));
    end

    %compute the theoretical allan variance for given model
    if (~isempty(t_axis))
        allan_std(in,:)=((qc*Tc)^2./t_axis.*(1-(Tc/2./t_axis).*(3-4*exp(-t_axis/Tc)+exp(-2*t_axis/Tc)))).^0.5;
    else
        allan_std=0;
    end
end

