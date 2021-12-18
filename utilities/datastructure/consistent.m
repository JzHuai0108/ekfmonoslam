function errstring = consistent(ds, type)

% CONSISTENT   Check ReBEL data structures for consistentency.
%
%   errstring = consistent(ds, type) checks a data structure DS of supposed TYPE for consistency.
%
%     See also
%     GSSM, GENINFDS
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


%===============================================================================================

if (nargin~=2)
  error('Not enough input arguments!');
end

if ~ischar(type)
  error('Type argument should be a string.');
end

% Assume that all is OK as default
errhead = ['  [ consistent(' inputname(1) ',''' type ''') ]   '];

errstring = '';


if ~isstruct(ds)
  errstring = [errhead inputname(1) ' is not a data structure.'];
  return
end


switch (type)

%===============================================================================================================================
%== General State Space Model Data Structure (GSSM)
%===============================================================================================================================
case 'gssm'

  %-- Check existence of required data structure fields
  mf = checkstructfields(ds,'type','tag','statedim','obsdim','paramdim','U1dim','U2dim','Vdim','Ndim','pNoise','oNoise','params','ffun','hfun','setparams');
  if ~isempty(mf), errstring = [errhead ' Data structure does not contain the following fields : ' mf]; return; end

  %-- Check type field
  if ~ischar(ds.type) | ~stringmatch(ds.type,'gssm')
      errstring = [errhead  ' Data structure is not of type ''' type ''' or the type field is not a string'];
      return;
  end


%===============================================================================================================================
%== Inference Data Structure (InferenceDS)
%===============================================================================================================================
case 'InferenceDS'

  %-- Check existence of required data structure fields
  mf = checkstructfields(ds,'type','inftype','statedim','obsdim','U1dim','U2dim','Vdim','Ndim', ...
                            'ffun','hfun','innovation','stateAngleCompIdxVec','obsAngleCompIdxVec');
  if ~isempty(mf), errstring = [errhead ' Data structure does not contain the following fields : ' mf]; return; end

  %-- Check type field
  if (~ischar(ds.type) | ~stringmatch(ds.type,'InferenceDS'))
    errstring = [errhead  ' Data structure is not of type ''' type ''' or the type field is not a string.'];
    return;
  end

  %-- Check inference type field
  if ~stringmatch(ds.inftype,{'state','parameter','joint'})
    errstring = [errhead  ' Inference type ''' type ''' unknown.'];
    return;
  end


%===============================================================================================================================
%== Noise Source Data Structure (NoiseDS)
%===============================================================================================================================
case 'NoiseDS'

  %-- Check existence of required data structure fields
  mf = checkstructfields(ds,'type','ns_type','tag','dim','sample','likelihood');
  if ~isempty(mf), errstring = [errhead ' Data structure does not contain the following fields : ' mf]; return; end

  %-- Check type field
  if (~ischar(ds.type) | ~stringmatch(ds.type,'NoiseDS'))
    errstring = [errhead  ' Data structure is not of type ''' type ''' or the type field is not a string.'];
    return;
  end


%===============================================================================================================================
otherwise
  errstring = [errhead 'Data structure type ' type ' not supported.'];
end

return
