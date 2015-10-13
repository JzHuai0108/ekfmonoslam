%sculling (essentially sculling+rotation) correction for each minor
%interval
%gyro & acc must be increment type (3,length) outputs. alg is the desired
%sculling algorithm
function [vinc, ainc, corr]=sculling_minor(gyro, acc, alg)
inlen=size(acc,2);
switch (alg)
    case (0)
        %ignagni(1998):algorithm 1 -> (1996):algo2
        outlen=floor((inlen-3)/2)+1;
        vinc=zeros(3,outlen);
        ainc=zeros(3,outlen);
        corr=zeros(3,outlen);

        ind=1;
        for i=3:2:inlen
            vinc(:,ind)=sum(acc(:,i-1:i),2);
            ainc(:,ind)=sum(gyro(:,i-1:i),2);
            
            rot=0.5*cross(ainc(:,ind), vinc(:,ind));
            vr_a=cross(((-1/30)*gyro(:,i-2)+(11/15)*gyro(:,i-1)), acc(:,i));
            vr_b=cross(((-1/30)*acc(:,i-2)+(11/15)*acc(:,i-1)), gyro(:,i));
            corr(:,ind)=rot+vr_a+vr_b;
            ind=ind+1;
        end
    case (1)
        %ignagni(1998):algorithm 2 -> (1996):algo3
        outlen=floor((inlen-3)/3)+1;
        vinc=zeros(3,outlen);
        ainc=zeros(3,outlen);
        corr=zeros(3,outlen);

        ind=1;
        for i=3:3:inlen
            vinc(:,ind)=sum(acc(:,i-2:i),2);
            ainc(:,ind)=sum(gyro(:,i-2:i),2);
            
            rot=0.5*cross(ainc(:,ind), vinc(:,ind));
            vr_a=cross(((9/20)*gyro(:,i-2)+(27/20)*gyro(:,i-1)), acc(:,i));
            vr_b=cross(((9/20)*acc(:,i-2)+(27/20)*acc(:,i-1)), gyro(:,i));
            corr(:,ind)=rot+vr_a+vr_b;
            ind=ind+1;
        end
    case(2)
        %ignagni(1998):algorithm 3 -> (1996):algo5
        outlen=floor((inlen-4)/3)+1;
        vinc=zeros(3,outlen);
        ainc=zeros(3,outlen);
        corr=zeros(3,outlen);

        ind=1;
        for i=4:3:inlen
            vinc(:,ind)=sum(acc(:,i-2:i),2);
            ainc(:,ind)=sum(gyro(:,i-2:i),2);
            
            rot=0.5*cross(ainc(:,ind), vinc(:,ind));
            vr_a=cross(((3/280)*gyro(:,i-3)+(57/140)*gyro(:,i-2)+(393/280)*gyro(:,i-1)), acc(:,i));
            vr_b=cross(((3/280)*acc(:,i-3)+(57/140)*acc(:,i-2)+(393/280)*acc(:,i-1)), gyro(:,i));
            corr(:,ind)=rot+vr_a+vr_b;
            ind=ind+1;
        end
    case(3)
        %ignagni(1998):algorithm 4 -> (1996):algo6
        outlen=floor((inlen-5)/3)+1;
        vinc=zeros(3,outlen);
        ainc=zeros(3,outlen);
        corr=zeros(3,outlen);

        ind=1;
        for i=5:3:inlen
            vinc(:,ind)=sum(acc(:,i-2:i),2);
            ainc(:,ind)=sum(gyro(:,i-2:i),2);
            
            rot=0.5*cross(ainc(:,ind), vinc(:,ind));
            vr_a=cross(((-1/420)*gyro(:,i-4)+(1/40)*gyro(:,i-3)+(157/420)*gyro(:,i-2)+(1207/840)*gyro(:,i-1)), acc(:,i));
            vr_b=cross(((-1/420)*acc(:,i-4)+(1/40)*acc(:,i-3)+(157/420)*acc(:,i-2)+(1207/840)*acc(:,i-1)), gyro(:,i));
            corr(:,ind)=rot+vr_a+vr_b;
            ind=ind+1;
        end
    case(4)
        %ignagni(1998):algorithm 5 -> (1996):algo7
        outlen=floor((inlen-4)/4)+1;
        vinc=zeros(3,outlen);
        ainc=zeros(3,outlen);
        corr=zeros(3,outlen);

        ind=1;
        for i=4:4:inlen
            vinc(:,ind)=sum(acc(:,i-3:i),2);
            ainc(:,ind)=sum(gyro(:,i-3:i),2);
            
            rot=0.5*cross(ainc(:,ind), vinc(:,ind));
            vr_a=cross(((54/105)*gyro(:,i-3)+(92/105)*gyro(:,i-2)+(214/105)*gyro(:,i-1)), acc(:,i));
            vr_b=cross(((54/105)*acc(:,i-3)+(92/105)*acc(:,i-2)+(214/105)*acc(:,i-1)), gyro(:,i));
            corr(:,ind)=rot+vr_a+vr_b;
            ind=ind+1;
        end
    case(5) %%quadratic fit to both gyro and accel outputs (similar to the coning algorithm in ignagni(1990):algoD
        outlen=floor((inlen-3)/3)+1;
        vinc=zeros(3,outlen);
        ainc=zeros(3,outlen);
        corr=zeros(3,outlen);
        
        
        ind=1;
        for i=3:3:inlen
            vinc(:,ind)=sum(acc(:,i-2:i),2);
            ainc(:,ind)=sum(gyro(:,i-2:i),2);
            
            %compute the polynomial coefficients
            b=(11/6)*gyro(:,i-2)+(-7/6)*gyro(:,i-1)+(1/3)*gyro(:,i);
            c=(-2)*gyro(:,i-2)+(3)*gyro(:,i-1)+(-1)*gyro(:,i);
            d=(1/2)*gyro(:,i-2)+(-1)*gyro(:,i-1)+(1/2)*gyro(:,i);
            
            e=(11/6)*acc(:,i-2)+(-7/6)*acc(:,i-1)+(1/3)*acc(:,i);
            f=(-2)*acc(:,i-2)+(3)*acc(:,i-1)+(-1)*acc(:,i);
            g=(1/2)*acc(:,i-2)+(-1)*acc(:,i-1)+(1/2)*acc(:,i);
            
            %Compute sculling(+rotation)
            corr(:,ind)=(3^2/2)*cross(b,e)+(3^3/3)*(cross(b,f)+cross(c,e)/2)+(3^4/4)*(cross(b,g)+cross(c,f)/2+cross(d,e)/3)+...
                (3^5/5)*(cross(c,g)/2+cross(d,f)/3)+(3^6/6)*cross(d,g)/3;
            
            ind=ind+1;
        end
        
    otherwise
        disp('undefined algorithm');
end



