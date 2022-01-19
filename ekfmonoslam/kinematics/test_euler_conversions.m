function test_euler_conversions()
eul = rand(3, 1);
R = euler2dcm_v000(eul);
eul2 = dcm2euler_v000(R);
eul3 = rotro2eu('xyz', R);
assert(norm(eul2 - eul3) < 3e-7);

eul4 = Cbn2att(R);
assert(norm(eul3 - eul4) < 3e-7);
end
