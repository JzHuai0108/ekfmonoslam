function Kss=cp_SSgain(A,Q,H,R,P, mtd)
num_st=size(A,1);
    
if (mtd==0)     %method to use when dare fails to reach a solution
    % [Beq, yy]=chol(Q,'lower');
    Beq=Q^0.5;
   
    %basis for observable space
    Ob=obsv(A,H);
    [u,s,v]=svd(Ob'*Ob);
    num_null_obs=length(find(diag(s)<1e-6));
    num_range_obs=num_st-num_null_obs;
    base_null_obs=v(:,(num_range_obs+1):num_st);
    base_range_obs=u(:,1:num_range_obs);
    Wo=[base_null_obs base_range_obs];   %x=Wz
    
    %basis for controllable and uncontrollable spaces
    Co=ctrb(A,Beq);
    [u,s,v]=svd(Co*Co');
    num_null_ctr=length(find(diag(s)<1e-12));
    num_range_ctr=size(A,1)-num_null_ctr;
    base_null_ctr=v(:,(num_range_ctr+1):size(A,1));
    base_range_ctr=u(:,1:num_range_ctr);
    Wc=[base_null_ctr base_range_ctr];   %x=Wz
    
    %%%compute the intersection of observable space with cont/uncont spaces
    %obs+ctr space (min real)
    if (num_null_ctr>0)
        mx_a=[base_range_obs base_range_ctr];
        [mx_b,pvt] = rref(mx_a);
        ind=setdiff(1:(num_range_obs+num_range_ctr),pvt);
        if (~isempty(ind))
            Woc=base_range_obs*mx_b(1:num_range_obs,ind);
            Woc*diag(diag(Woc'*Woc).^-0.5); %normalize
            num_base_oc=size(Woc,2);
        else
            Woc=[];
            num_base_oc=0;
        end
        %obs+unctr space
        mx_a=[base_range_obs base_null_ctr];
        [mx_b,pvt] = rref(mx_a);
        ind=setdiff(1:(num_range_obs+num_null_ctr),pvt);
        if (~isempty(ind))
            Wou=base_range_obs*mx_b(1:num_range_obs,ind);
            Wou=diag(diag(Wou'*Wou).^-0.5); %normalize
            num_base_ou=size(Wou,2);
        else
            Wou=[];
            num_base_ou=0;
        end
        
        %%%overall transformation matrix
        W=[base_null_obs Wou Woc];
    else    %no uncontrollable subspace
        W=[base_null_obs base_range_obs];
        num_base_ou=0;
        num_base_oc=base_range_obs;
    end
    
    %transform the system
    W_inv=inv(W);
    A_t=W_inv*A*W;
    Beq_t=W_inv*Beq;
    Qeq_t=Beq_t*Beq_t';
    H_t=H*W;
    
    %partitions
    in1=num_null_obs;
    in2=num_null_obs+num_base_ou;
    in3=num_null_obs+num_base_ou+num_base_oc;
    
    Au=A_t(1:in1,1:in1);   %unobservable part
    Au_cross=A_t(1:in1,in1+1:end);   %cross effects of observable part on unobservable part
    Aou=A_t(in1+1:in2,in1+1:in2);   %observable but unctrollable part
    Aoc=A_t(in2+1:end,in2+1:end);   %observable and controllable
    Aoc_cross=A_t(in2+1:end,in1+1:in2); %cross effects of ou part on oc part (not used)
    Ao=A_t(in1+1:end,in1+1:end);   %observable part
    
    Qu=Qeq_t(1:in1,1:in1);
    Qu_cross=Qeq_t(1:in1,in1+1:end);    %cross disturbance between o and uno parts
    Qoc=Qeq_t(in2+1:end,in2+1:end);
    Qo=Qeq_t(in1+1:end,in1+1:end);
    
    Hou=H_t(:,in1+1:in2);
    Hoc=H_t(:,in2+1:end);
    Ho=H_t(:,in1+1:end);    %observable
    
    %steady state solution for ubsorvable and controllable parts (including
    %uncont. parts except the unit circle)
    [Poc,clp,Koc]=dare(Aoc',Hoc', Qoc, R);
    Koc=Poc*Hoc'*inv(Hoc*Poc*Hoc'+R);
    
    %Add SS gain values for the unit circle uncontrollable parts but observable
    %parts
    Ko=[zeros(num_base_ou,size(Koc,2));Koc];
    Po=diag_mat(Poc,1,zeros(num_base_ou));
    
    % %%Steady state solution for observable part
    % [Po,clp,Ko]=dare(A3',H2', Q3, R);
    % Ko=Po*H2'*inv(H2*Po*H2'+R);
    
    %%%%Compute the cross gain
    mx_a=Au;
    mx_b=(eye(num_range_obs)-Ho'*inv(Ho*Po*Ho'+R)*Ho*Po)*Ao';
    mx_c=Au_cross*(Po-Ko*Ho*Po)'*Ao'+Qu_cross;
    %%check the multiplication of eigenvalues
    if (isempty(mx_a)) %no unobservable space
        eig_mul=0;
    else
        lef=eig(mx_a);
        rig=eig(mx_b);
        eig_mul=lef*rig';
    end
    
    if (isempty(find(eig_mul==1,1,'first')))
        Pu_cross = dlyap(mx_a,mx_b,mx_c);
        Ku_cross=Pu_cross*Ho'*inv(Ho*Po*Ho'+R);
        %%overall gain
        Kss=W*[Ku_cross;Ko];
    else %in this case stein equation has no unique solution (solution depends initial P value)
        %Perform standard KF iteration to find a K
        disp('comp_ss_gain: NO SS!!');
        K_pre=zeros(num_st,size(H,1));
        K=P*H'*inv(H*P*H'+R);
        while (max(abs(K(:)-K_pre(:))>1e-7))
            K_pre=K;
            P=(eye(size(A))-K*H)*P;
            P=A*P*A'+Beq*Beq';
            K=P*H'*inv(H*P*H'+R);
        end
        Kss=K;
    end
elseif (mtd==1) %iterations
    K_pre=zeros(num_st,size(H,1));
    K=P*H'*inv(H*P*H'+R);
    in=0;
    while (max(abs(K(:)-K_pre(:))>1e-7))
        K_pre=K;
        P=(eye(size(A))-K*H)*P;
        P=A*P*A'+Q;
        K=P*H'*inv(H*P*H'+R);
        in=in+1;
    end
    disp(['cpSSgain 145: SS iteration =' num2str(in)]);
    Kss=K;
end