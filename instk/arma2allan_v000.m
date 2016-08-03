%assumption: number of zeros must be strictly less than number of poles.
%(a0)x(n+m)+(a1)x(n+m-1)...+(am)x(n)=(b0)w(n+m)+(b1)w(n+m-1)...+(bk)a(n+m-k
%) 
function allan_th=arma2allan(m_A, m_B, fs,t_axis, allan)

[nmod na]=size(m_A);
nb=size(m_B,2);
allan_th=zeros(nmod,length(t_axis));

%discrete indeces
n_axis=round(t_axis*fs);

%normalize the model so that m_A(1)=1
m_A=m_A./(m_A(:,1)*ones(1, size(m_A,2)));
m_B=m_B./(m_A(:,1)*ones(1, size(m_B,2)));


for in=1:nmod
    %stip the excessive zeros for each model
    A=m_A(in,:);
    B=m_B(in,:);
    sr_a=find(fliplr(A)~=0,1);
    if (~isempty(sr_a))
        A=A(1,1:end-sr_a+1);
    end
    sr_a=find(fliplr(B)~=0,1);
    if (~isempty(sr_a))
        B=B(1,1:end-sr_a+1);
    end
    
    if (size(A,2)==1) %WN or MA
        if (size(B,2)==1) %WN
            allan_th(in,:)=m_B(in,1)./((n_axis).^0.5);
        else %MA
            disp('error - To be completed later');
			%%%%TOOOO LAZZZZZYYYY
        end
    elseif ((size(A,2)==2) && A(2)==(-1)) %RW
        allan_th(in,:)=((2*(n_axis.^2)+1)./(6*n_axis)).^0.5*B(1);
    else %strictly proper ARMA
        %convert model into correlation exponentials
        exp_param=arma2cov_v000(B, A);
        allan_th(in,:)=exp2allan(exp_param,n_axis);
    end
end

% if (~isempty(allan))
%     figure;
%     loglog(t_axis,allan);grid;
%     hold on;
%     for in=1:nmod
%         loglog(t_axis,allan_th(in,:),'k');
%     end
%     loglog(t_axis,sum(allan_th.^2,1).^0.5,'r');
% end

if (~isempty(allan))
    figure;
    loglog(t_axis,allan,'linewidth',2);
    hold on;
    loglog(t_axis,sum(allan_th.^2,1).^0.5,'r','linewidth',2);   %%just for legend correction
    for in=1:nmod
        loglog(t_axis,allan_th(in,:),'k','linewidth',2);
    end
    loglog(t_axis,sum(allan_th.^2,1).^0.5,'r','linewidth',2);
end
xlabel('Cluster Length (sec)');
%ylabel('Allan Std (rad/sec)');
grid;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function allan_res=exp2allan(ep,n_axis)
nexp=size(ep,1);
allan_res=zeros(1,length(n_axis));
for in=1:length(n_axis)
    n=n_axis(in);
    var_cluster=sum(ep(:,2))/n;
    cor_cluster=0;
    for ii=1:nexp
        c=ep(ii,2);
        v=ep(ii,1);
        sr_a=2/n*c*v*(n-1-n*v+v^n)/n/(1-v)^2;
        var_cluster=var_cluster+sr_a;

        sr_a=c*(v^n)/n;
        sr_b=c/(n^2)*v*(1-n*(v^(n-1))+(n-1)*(v^n))/(1-v)^2;
        sr_c=c/(n^2)*v^(2*n-1)*(1-n*(v^(1-n))+(n-1)*(v^(-n)))/(1-v^(-1))^2;
        
        sr_d=sr_a+sr_b+sr_c;
        if (isnan(sr_d))    %this method is very problematic in terms of numerical precison for high "n" values
            cor_cluster=cor_cluster+0;
        else
            cor_cluster=cor_cluster+sr_d;
        end
    end

    allan_res(in)=(var_cluster-cor_cluster)^0.5;
end
    
            