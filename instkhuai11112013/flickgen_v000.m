%use nstage number of poles and zeros to simulate flicker noise (h/2f) between
%fmin and fmax. To obtain acceptable approximation nstage must be equal or bigger
%than log10(fmax/fmin). The bigger the nstage, the better the approximation
%for the middle part of the frequency range (however the approximation will
%be worse for the regions around fmin and fmax.)

function [gain, pls, zrs]=flickgen(h, fmin, fmax, nstage)
wmin=fmin*2*pi;
wmax=fmax*2*pi;

if isempty(nstage)
    nstage=floor(log10(fmax/fmin));
end
wstep=(wmax/wmin)^(1/(4*nstage));
wmid=sqrt(wmax*wmin);
pls=zeros(nstage,1);
zrs=zeros(nstage,1);

%pol/zero locations (placed evenly in the logaritmic scale) and the gain
magmid=1;
for i=1:nstage
    pls(i)=-(wmin*wstep^(4*i-3));
    zrs(i)=-(wmin*wstep^(4*i-1));
    magmid=magmid*(wmid^2+zrs(i)^2)/(wmid^2+pls(i)^2);
end
gain=sqrt(h*pi/wmid/magmid);

return;

%%%%%%%%%%%% Verbose %%%%%%%
%%Compare the psd of approximation to flicker noise
w=logspace(log10(wmin/10), log10(wmax*10),1000);
mag=1;
for i=1:nstage
    mag=mag.*(w.^2+zrs(i).^2)./(w.^2+pls(i).^2);
end
mag=gain^2*mag;

loglog(w,mag);
hold on
loglog(wmin:wmax, h*pi./(wmin:wmax),'r');
grid;

%%Compare the allan variance of approximation to flicker noise
tau=logspace(log10(1/fmax/2), log10(1/fmin*2), 50);
allan=zeros(length(tau),1);
dw=diff(w);
for i=1:length(tau)
    allankernel=sinc(w*tau(i)/2/pi).^2.*sin(w*tau(i)/2).^2;
    spect=allankernel.*mag*2/pi;
    allan(i)=sum(dw.*spect(1:end-1));
end

figure;
loglog(tau, allan);
hold on
loglog(tau, h*log(2)*2,'r');

