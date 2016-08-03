%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fclose all;
clear;
PATH='C:\Work\Benisil';
nst=15;

%this structure is coded for the filtered estimates (not the predicted).
%There are infinetely many different ways of implementation!!!! This is
%just an example.

%fixed interval smoother operations
F_RES=fopen([PATH '\res1.bin'],'rb');   %Filter results
F_PROP=fopen([PATH '\smo_prop.bin'],'rb');  %STM and covariances
F_UPDT=fopen([PATH '\smo_updt.bin'],'rb');  %innovations for each update
F_SMO=fopen([PATH '\res_bf.bin'],'wb');

status=fseek(F_RES, 0, 'eof');
status=fseek(F_PROP, 0, 'eof');
status=fseek(F_UPDT, 0, 'eof');

sz_res=(2*nst)+1;
sz_prop=(2*nst^2)+1;
sz_updt=(2*nst^2+nst)+1;

dx_accum=zeros(nst,1);
dP_accum=zeros(nst);

%%read the last results
fseek(F_RES, -sz_res*8, 'cof');
fseek(F_PROP, -sz_prop*8, 'cof');
fseek(F_UPDT, -sz_updt*8, 'cof');
f_res=fread(F_RES,sz_res,'double');
f_prop=fread(F_PROP,sz_prop,'double');
f_updt=fread(F_UPDT,sz_updt,'double');
pr_coun=f_prop(1);
up_coun=f_updt(1);

while (pr_coun~=-1)
    P_filt=reshape(f_prop(2:1+nst^2),nst,nst);
    STM_prop=reshape(f_prop((2+nst^2):end),nst,nst);
    filt_err=f_res(2:nst+1);
    
    %Smoother estimates (computed based on filtered results)
    smo_err=filt_err-P_filt*dx_accum;    %smo_err is computed via updated-filter states
    P_smo=(eye(nst)-P_filt*dP_accum)*P_filt;  %smo covariance
    
    %write the results to the file
    fwrite(F_SMO,[pr_coun;smo_err;diag(P_smo).^0.5],'double');
    
    %if there are new observations, update the smoother states
    if (up_coun==pr_coun)  %in case there is a update      
        dx_inp=f_updt(2:nst+1);
        dP_inp=reshape(f_updt(2+nst:1+nst+nst^2),nst,nst);

        STM_upd=reshape(f_updt(2+nst+nst^2:end),nst,nst);
        dx_accum=STM_upd'*dx_accum+dx_inp;   %next value of dx_accum 
        dP_accum=STM_upd'*dP_accum*STM_upd+dP_inp;
        
        %read the next value of f_updt
        status=fseek(F_UPDT, -sz_updt*8*2, 'cof');
        if status==-1
            up_coun=-1;
        else
            f_updt=fread(F_UPDT,sz_updt,'double');
            up_coun=f_updt(1);
        end
    end
    
    %Propagate the smoothed states to the next output time
    dx_accum=STM_prop'*dx_accum;   %next value of dx_accum 
    dP_accum=STM_prop'*dP_accum*STM_prop;
       
    %read the values at the next propagation cycle.
    status=fseek(F_RES, -sz_res*8*2, 'cof');
    status=fseek(F_PROP, -sz_prop*8*2, 'cof');
    if (status==-1)
        pr_coun=-1;
    else
        f_res=fread(F_RES,sz_res,'double');
        f_prop=fread(F_PROP,sz_prop,'double');
        pr_coun=f_prop(1);
    end
end

fclose all

%plot the results
figure;
res1=readbin_v000([PATH '\res1.bin'],2*nst+1);
ress=readbin_v000([PATH '\res_bf.bin'],2*nst+1);
ress=fliplr(ress);
plot(res1(1,:)/100,abs(res1(2,:))*6e6)
hold on
plot(ress(1,:)/100,abs(ress(2,:))*6e6,'r')
plot(res1(1,:)/100,abs(res1(nst+2,:))*6e6,'--')
plot(ress(1,:)/100,abs(ress(nst+2,:))*6e6,'--r')
title('lat. error');
ylabel('m');
xlabel('sec');
legend('filter error', 'bf (smoother) error', 'filter std','bf cov');
grid;
return;

figure;
yy=5+1;
plot(res1(1,:),abs(res1(yy,:)))
hold on
plot(ress(1,:),abs(ress(yy,:)),'r')
plot(res1(1,:),abs(res1(nst+yy,:)),'--')
plot(ress(1,:),abs(ress(nst+yy,:)),'--r')
xlabel('sec');
legend('filter', 'bf');
grid;
return;