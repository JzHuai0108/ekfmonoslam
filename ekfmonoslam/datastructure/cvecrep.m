function m = cvecrep(v,c)

% CVECREP  Column vector replicate
%
%   M = cvecrep(V, C) Replicates a Nx1 dimensional column vector V, C times to generate a
%   NxC dimensional matrix M.
%
%   See also
%   RVECREP, REPMAT
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

if isempty(v)

  m = zeros(0,c);

else

  m = v(:,ones(c,1));

end
