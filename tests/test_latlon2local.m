% check latlon2local

origin = [30.5298813668333; 114.350355321667; 20.8869];
lla = [30.5298266666667; 114.350367833333; 20.8];
lla = [30.5299376666667; 114.35528; 17.6];
% matlab functions
[e1, n1, u1] = geodetic2enu(lla(1), lla(2), lla(3), origin(1),...
    origin(2), origin(3), wgs84Ellipsoid, 'degrees');
[e2, n2, u2] = latlon2local(lla(1), lla(2), lla(3), origin);

assert(abs(e1 - e2) < 1e-8);
assert(abs(n1 - n2) < 1e-8);
assert(abs(u1 - u2) < 1e-8);

% our own implementation
Ce2n0=llh2dcm_v000(origin(1:2) * pi/ 180,[0;1]);
qs02e=rotro2qr(Ce2n0');
origin_xyz = ecef2geo_v000([origin(1:2) * pi / 180; origin(3)],1);
xyz = ecef2geo_v000([lla(1:2) * pi / 180; lla(3)],1);
ned = quatrot_v000(qs02e, xyz - origin_xyz ,1)';

assert(abs(e1 - ned(2)) < 1e-8);
assert(abs(n1 - ned(1)) < 1e-8);
assert(abs(u1 + ned(3)) < 1e-8);