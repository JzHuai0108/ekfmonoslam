%RPY2TR	Roll/pitch/yaw to homogenous transform
%
%	RPY2TR([R P Y])
%	RPY2TR(R,P,Y) returns a homogeneous tranformation for the specified
%	roll/pitch/yaw angles.  These correspond to rotations about the
%	Z, X, Y axes respectively.
%
%	See also TR2RPY, EUL2TR

%	Copright (C) Peter Corke 1993
function r = rpy2tr(roll, pitch, yaw)
        if length(roll) == 3,
                r = rotz(roll(1)) * roty(roll(2)) * rotx(roll(3));
        else
                r = rotz(roll) * roty(pitch) * rotx(yaw);
        end
