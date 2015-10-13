%% Input is the rate output and frequency in Hertz
%rate_data=3*ndat (col mat), max_step in minutes
function [t_axis, allan, Nrw]=cp_allan_var(rate_data, fs, show_fig, max_step)
[naxe ndat]=size(rate_data);
clus_max=floor((ndat-1)/2);

if (isempty(max_step))  %not a compulsory input
    max_step=clus_max;
else
    max_step=max_step*60*fs; %convert from minutes to number of samples
    if (max_step>clus_max)
        max_step=clus_max;
    end
end

%use log2 interval upto max_step
sr_a=floor(log2(max_step));
clus_size=2.^(0:sr_a);

%use linear interval for the rest
clus_size=[clus_size max_step:max_step:clus_max];

%total number of clusters
clus_num=length(clus_size);

%time axis corresponding to the given cluster lengths
t_axis=clus_size/fs;

%Compute allan variance
int_data=cumsum(rate_data,2);
allan=zeros(naxe, clus_num);
for i=1:clus_num
    n=clus_size(i);      %Length of each cluster
    vr_c=int_data(:,1+2*n:ndat); 
    vr_b=int_data(:,1+n:ndat-n);
    vr_a=int_data(:,1:ndat-2*n);
    allan(:,i) = sqrt(mean((vr_c-2*vr_b+vr_a).^2,2)) /sqrt(2) /n;
end

n=round(log2(fs))+1;
Nrw=allan(:,n).^2*t_axis(n); %rw power (=cont. time psd level)

if (show_fig~=0)
    ccodes=['b','r','k','m','g','c','y'];
    figure(show_fig);
    for in=1:naxe
        arw_line=(Nrw(in)./t_axis).^0.5;
        loglog(t_axis,allan(in,:),['-d' ccodes(in)]);
        hold on;
        loglog(t_axis,arw_line,'.r');
        %legend('Allan Std','ARW Line');
    end
    xlabel('Cluster Length (Sec)');
    grid;
end