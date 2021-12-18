function m = rvecrep(v,c)

% RVECREP  Row vector replicate
%
%   M = rvecrep(V, C) Replicates a 1xN dimensional row vector V, C times to generate a
%   CxN dimensional matrix M.
%
%   See also
%   CVECREP, REPMAT

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

    m = zeros(c,0);

else

    m = v(ones(1,c),:);

end
