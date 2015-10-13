function res=skew(inp_vr)
res=[0 			-inp_vr(3)	inp_vr(2); ...
     inp_vr(3)	0			-inp_vr(1);...
	 -inp_vr(2)	inp_vr(1)	0];
	 

%res=zeros(3,3);
%res(1,2)=-inp_vr(3);
%res(1,3)=inp_vr(2);
%res(2,1)=inp_vr(3);
%res(2,3)=-inp_vr(1);
%res(3,1)=-inp_vr(2);
%res(3,2)=inp_vr(1);

