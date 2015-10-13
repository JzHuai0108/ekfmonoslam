%data   :3*length increment type equally spaced sensor output.
%alg    :Type of the coning algorithm to be used   
function [inc corr]=coning_minor(data, alg)
inlen=size(data,2);
switch (alg)
    case(0)
        %ignagni(1990):Algorithm A (no correction - constant approximation)     
        inc=data;
        corr=zeros(size(data));
    case(1)
        %ignagni(1990):Algorithm D (quadratic approximation to 3 points)    
        outlen=floor((inlen-3)/3)+1;
        inc=zeros(3,outlen);
        corr=zeros(3,outlen);
        
        ind=1;
        for i=3:3:inlen
            inc(:,ind)=sum(data(:,i-2:i),2);
            corr(:,ind)=(33/80)*cross(data(:,i-2), data(:,i))+(57/80)*cross(data(:,i-1),data(:,i)-data(:,i-2));
            ind=ind+1;
        end
    case(2)
        %ignagni(1990):Algorithm E
        outlen=floor((inlen-3)/3)+1;
        inc=zeros(3,outlen);
        corr=zeros(3,outlen);
        
        ind=1;
        for i=3:3:inlen
            inc(:,ind)=sum(data(:,i-2:i),2);
            corr(:,ind)=(9/20)*cross(data(:,i-2), data(:,i))+(27/40)*cross(data(:,i-1),data(:,i)-data(:,i-2));
            ind=ind+1;
        end
    case(3)
        %ignagni(1996):Algorithm 1
        outlen=floor((inlen-4)/2)+1;
        inc=zeros(3,outlen);
        corr=zeros(3,outlen);
        
        ind=1;
        for i=4:2:inlen
            vr_a=sum(data(:,i-1:i),2);
            vr_b=sum(data(:,i-3:i-2),2);
            inc(:,ind)=var_a;
            corr(:,ind)=(32/45)*cross(data(:,i-1), data(:,i))+(-1/180)*cross(vr_b,vr_a);
            ind=ind+1;
        end
    case(4)
        %ignagni(1996):Algorithm 2
        outlen=floor((inlen-3)/2)+1;
        inc=zeros(3,outlen);
        corr=zeros(3,outlen);

        ind=1;
        for i=3:2:inlen
            inc(:,ind)=sum(data(:,i-1:i),2);
            corr(:,ind)=cross(((-1/30)*data(:,i-2)+(11/15)*data(:,i-1)), data(:,i));
            ind=ind+1;
        end
    case (5)
        %ignagni(1996):Algorithm 3 - Ignagni (1990):Algorithm F
        outlen=floor((inlen-3)/3)+1;
        inc=zeros(3,outlen);
        corr=zeros(3,outlen);

        ind=1;
        for i=3:3:inlen
            inc(:,ind)=sum(data(:,i-2:i),2);
            corr(:,ind)=cross(((9/20)*data(:,i-2)+(27/20)*data(:,i-1)), data(:,i));
            ind=ind+1;
        end
    case (6)
        %ignagni(1996):Algorithm 4
        outlen=floor((inlen-6)/3)+1;
        inc=zeros(3,outlen);
        corr=zeros(3,outlen);

        ind=1;
        for i=6:3:inlen
            vr_a=sum(data(:,i-2:i),2);
            vr_b=sum(data(:,i-5:i-3),2);
            inc(:,ind)=vr_a;
            corr(:,ind)=cross(((243/560)*data(:,i-2)+(1539/1120)*data(:,i-1)), data(:,i))+(1/3360)*cross(vr_b,vr_a);
            ind=ind+1;
        end
    case (7)
        %ignagni(1996):Algorithm 5
        outlen=floor((inlen-4)/3)+1;
        inc=zeros(3,outlen);
        corr=zeros(3,outlen);

        ind=1;
        for i=4:3:inlen
            inc(:,ind)=sum(data(:,i-2:i),2);
            corr(:,ind)=cross(((3/280)*data(:,i-3)+(57/140)*data(:,i-2)+(393/280)*data(:,i-1)), data(:,i));
            ind=ind+1;
        end
    case (8)
        %ignagni(1996):Algorithm 6
        outlen=floor((inlen-5)/3)+1;
        inc=zeros(3,outlen);
        corr=zeros(3,outlen);

        ind=1;
        for i=5:3:inlen
            inc(:,ind)=sum(data(:,i-2:i),2);
            corr(:,ind)=cross(((-1/420)*data(:,i-4)+(1/40)*data(:,i-3)+(157/420)*data(:,i-2)+(1207/840)*data(:,i-1)), data(:,i));
            ind=ind+1;
        end   
    case (9)
        %ignagni(1996):Algorithm 7
        outlen=floor((inlen-4)/4)+1;
        inc=zeros(3,outlen);
        corr=zeros(3,outlen);

        ind=1;
        for i=4:4:inlen
            inc(:,ind)=sum(data(:,i-3:i),2);
            corr(:,ind)=cross(((54/105)*data(:,i-3)+(92/105)*data(:,i-2)+(214/105)*data(:,i-1)), data(:,i));
            ind=ind+1;
        end
    otherwise
        disp('Undefined Algorithm (Perhaps, someone else adds it later)');
end