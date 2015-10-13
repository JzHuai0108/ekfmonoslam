function [pos_e, vel_e, Cbn_e, rstat]=AddIniErr(pos, vel, Cbn, errdef, rstat)

%randomize
if (rstat(1)~=0)
    randn('state',rstat)
else
    rstat=randn('state');
end

if (~isempty(errdef))
    pos_e=pos+errdef.pos_sP*randn(3,1);
    vel_e=vel+errdef.vel_sP*randn(3,1);
    att_e=errdef.att_sP*randn(3,1);

    Cerr=euler2dcm_v000(att_e);
    Cbn_e=Cerr'*Cbn;
else
    pos_e=pos;
    vel_e=vel;
    Cbn_e=Cbn;
end