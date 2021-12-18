%optimal solutions based on Wahba's problem.
%Cbn_est=(I-S(e))Cbn, e=E[b_errors]
%body=Cnb*nav+v;

function [Cnb E err]=align_opt(body, nav, coef)
    if isempty(coef)
        coef=ones(1, size(body,2));
    end

    if (length(coef)==2 && (coef(1)==0 || coef(2)==0)) %use algebraic method
        Cnb=algebraic(body,nav,coef);
        
        %%Compute the sensitivity matrix using the numeric methods
        E=numeric_sensitivity(@algebraic, body, nav,coef);
        
        mx_a=body-Cnb*nav;
        err=trace(mx_a'*mx_a);
    else %use the svd method
        Cnb=markley(body,nav,coef);
        
%         %%Compute the sensitivity matrix using the numeric methods
%         E=numeric_sensitivity(@markley, Cnb*nav,nav,coef);
        
        %%Compute the approximate sensitivity using algebraic method
        E=markley_sensitivity(Cnb*nav, nav, coef);
        
        mx_a=(body-Cnb*nav);
        err=coef*diag(mx_a'*mx_a);
    end
end

function Cnb=algebraic(body, nav, coef)
    %Determine the vector with 0 weight
    in1=1;
    in2=2;
    if (coef(1)==0)
        in1=2;
        in2=1;
    end
    b1=body(:,in1)/norm(body(:,in1));
    b2=body(:,in2)/norm(body(:,in2));
    n1=nav(:,in1)/norm(nav(:,in1));
    n2=nav(:,in2)/norm(nav(:,in2));

    vr_a=cross(b1, b2)/norm(cross(b1, b2));
    U=[b1 vr_a cross(b1,vr_a)];

    vr_a=cross(n1, n2)/norm(cross(n1, n2));
    V=[n1 vr_a cross(n1,vr_a)];

    Cnb=U*V';
end

%SVD based optimal solution to wahba's problem
function Cnb=markley(body,nav,coef)
    B=((ones(3,1)*coef).*body)*nav';
    [u,s,v]=svd(B);
    d=det(u)*det(v);
    Cnb=u*diag([1,1,d])*v';
    
    %lost function error (not the error on the original cost)
    %s(3,3)=d*s(3,3);
    %err=1-trace(s);
end

%Quaternion based approximate solution to wahba's problem
function [Cnb err]=schuster(body,nav,coef)
    %%todo
    Cnb=0;
    err=0;
end

%%Approximate sensitivity
function E=markley_sensitivity(body, nav, coef)
    B=((ones(3,1)*coef).*body)*nav';
    [u,s,v]=svd(B);
    d1=det(u);
    d2=det(v);
    U=u*diag([1 1 d1]);
    V=v*diag([1 1 d2]);
    Dinv=diag([1/(s(2,2)+d1*d2*s(3,3)), 1/(s(1,1)+d1*d2*s(3,3)), 1/(s(1,1)+s(2,2))]);
    Cnb=U*V';
    
    E=zeros(3,3*size(body,2));
    for in=1:size(body,2)
        E(:,in*3-2:in*3)=-coef(in)*skew(V'*nav(:,in))*U';
    end
    
    E=U*Dinv*E;
    E=-Cnb'*E;  %%Change for the notation
end


%%Numeric sensitivty
function E=numeric_sensitivity(fhandle, body,nav,coef)
    %nominal val
    Cbn=fhandle(body,nav,coef)';
    E=zeros(3,3*size(body,2));
    eps=0.001;
    for in=1:size(body,2)
        for i=1:3
            body_d=body;
            body_d(i,in)=body_d(i,in)+eps;
            Cnb_d=fhandle(body_d,nav,coef);
            mx_a=Cbn*Cnb_d-eye(3);
            E(:,i+(3*in-3))=[-mx_a(2,3);mx_a(1,3);-mx_a(1,2)]/eps;
        end
    end
end
