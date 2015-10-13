%given a direction vector computes the euler angles (assuming roll=0)
function res=dir2euler(dirv)
nv=size(dirv,1);
res=zeros(nv,3);
for in=1:nv
    sr_a=norm(dirv(in,:));
    pitch = asin(-dirv(in,3)/sr_a);
    head = atan2(dirv(in,2),dirv(in,1));
    res(in,:)=[0 pitch head];
end