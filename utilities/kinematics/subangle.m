function C = subangle(A, B)

%  ADDANGLE   Addition function for 'angle space' sigma-points expressed in radians.
%             This needed to deal with the angular discontinuety at +- pi radians.
%
%             C = addangle(A,B)
%
%   INPUT
%           A and B  : angles expressed in radians
%   OUTPUT
%           C        : sum C=A+B such that  -pi < C < pi
%
%   Copyright (c) Oregon Health & Science University (2006)
%
%   This file is part of the ReBEL Toolkit. The ReBEL Toolkit is available free for
%   academic use only (see included license file) and can be obtained from
%   http://choosh.csee.ogi.edu/rebel/.  Businesses wishing to obtain a copy of the
%   software should contact rebel@csee.ogi.edu for commercial licensing information.
%
%   See LICENSE (which should be part of the main toolkit distribution) for more
%   detail.

%=============================================================================================

C = A - B;

twopi = 2*pi;

idx1 = C > pi;
idx2 = C < -pi;

C(idx1) = C(idx1) - twopi;
C(idx2) = C(idx2) + twopi;

