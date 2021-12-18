function NoiseDS = gennoiseds(ArgDS)

% GENNOISEDS    Generates a NoiseDS data structure describing a noise source.
%
%   NoiseDS = gennoiseds(ArgDS)
%
%   This function generates a ReBEL noise source which is encapsulated in the NoiseDS data structure.
%   All arguments to the function are passed via the ArgDS argument data structure which has the
%   following fields, depending on the required noise source type:
%
%   1) Gaussian Noise Source   :  Single Gaussian noise source
%
%      ArgDS.type    : (string)    'gaussian'
%           .cov_type : (string)    'full'       full covariance matrix (default value)
%                                  'diag'       single diagonal covariance
%                                  'sqrt'       square-root form (lower triangular Cholesky factor)
%                                  'sqrt-diag'  square-root single diagonal form
%           .tag     : (string)     ID tag      (default='')
%           .dim     : (scalar)     noise vector length
%           .mu      : (c-vector)   mean vector
%           .cov     : (matrix)     covariance matrix (should comply with cov_type)
%
%   2) Combination Gaussian Noise Source : Combination of N independent Gaussian sources. All Gaussian
%                                          sources must have the same cov_type.
%
%      ArgDS.type     : (string)      'combo-gaussian'
%           .tag      : (string)      ID tag
%           .dim      : (scalar)      noise vector length
%           .noiseSources : (cell array)  cell array of N Guassian noise sources (NoiseDS data structures)
%
%
%   3) Generic Combination Noise Source : Combination of N independent noise sources.
%
%      ArgDS.type     : (string)      'combo'
%           .tag      : (string)      ID tag
%           .dim      : (scalar)      noise vector length
%           .noiseSources : (cell array)  cell array of N noise sources (NoiseDS data structure)
%
%
%   4) Scalar Gamma noise source
%
%      ArgDS.type     : (string)      'gamma'
%           .tag      : (string)      ID tag
%           .dim      : (scalar)      1 (multivariate Gamma noise not yet supported)
%           .alpha    : (scalar)      alpha parameter
%           .beta     : (scalar)      beta parameter
%
%
%   5) Gaussian Mixture noise Source
%
%      ArgDS.type     : (string)   'gmm'
%           .cov_type  : (string)   'full'       full covariance matrix (default value)
%                                  'diag'       single diagonal covariance
%                                  'sqrt'       square-root form (lower triangular Cholesky factor)
%                                  'sqrt-diag'  square-root single diagonal form
%           .tag     : (string)     ID tag      (default='')
%           .dim     : (scalar)     noise vector dimension
%           .M       : (scalar)     number of Gaussian component densities
%           .mu      : (dim-by-M matrix)  noise mean vectors of M mixture components
%           .cov     : (dim-by-dim-by-M matrix)  noise covariance matrices of N mixture components (should comply with cov_type)
%           .weights : (1-by-N vector)    mixing weights (priors) of mixture components
%
%
%   The generated noise source data structure, NoiseDS, has the following required fields. Depending on the noise source
%   type, the data structure may also have other type dependent fields. See documentation for full details.
%
%     NoiseDS.type       : (string) data structure type : 'NoiseDS'
%            .ns_type    : (string) noise source type
%            .dim        : (scalar) noise source dimension
%            .sample     : (function handle) Function to generate N noise source samples
%            .likelihood : (function handle) Function to evaluate the likelihood of a given noise source sample
%            .update     : (function handle) <<optional>> Function to update the internal structure of noise source.
%                                            This is not required of all noise sources.
%
%   See also
%     CONSISTENT, GENSYSNOISEDS
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

%========================================================================================================

%--- ERROR CHECKING & TYPE INDEPENDENT STRUCTURE ASSIGNMENT ---------------------------------------------

if (nargin ~= 1), error(' [ gennoiseds ] Incorrect number of inputs'); end

if ~isstruct(ArgDS), error(' [ gennoiseds ] The input argument to this function must be an argument data structure.'); end

if ~(isfield(ArgDS,'type') & ischar(ArgDS.type))
  error(' [ gennoiseds ] The argument data structure must have a ''type'' field, specifying the desired noise source type. This field must be a string.');
else
  Noise.type = ArgDS.type;       % assign noise source type
end

if isfield(ArgDS,'tag')          % assign tag
  if ischar(ArgDS.tag)
    Noise.tag = ArgDS.tag;
  else
    error(' [ gennoiseds ] The ''tag'' field must be a string.');
  end
else
  Noise.tag = '';
end

if (isfield(ArgDS,'dim') & isnumeric(ArgDS.dim) & (length(ArgDS.dim(:))==1))   %-- Check for and assign noise source dimension
  Noise.dim  = ArgDS.dim;
else
  error(' [ gennoiseds ] Noise source dimension not specified or not a scalar.');
end


%--- BUILD NOISE SOURCE DEPENDENT ON TYPE ----------------------------------------------------

switch (Noise.type)

%=============================================================================================
case 'gamma'

    if (Noise.dim~=1)
        error(' [ gennoiseds::gamma ] Multivariate Gamma noise not supported yet.');
    end

    if isfield(ArgDS,'alpha') & isnumeric(ArgDS.alpha) & (size(ArgDS.alpha) == [1 1])
        Noise.alpha = ArgDS.alpha;
    else
        error(' [ gennoiseds::gamma ] Alpha parameter must be a scalar');
    end

    if isfield(ArgDS,'beta') & isnumeric(ArgDS.beta) & (size(ArgDS.beta) == [1 1])
        Noise.beta = ArgDS.beta;
    else
        error(' [ gennoiseds::gamma ] Beta parameter must be a scalar');
    end

    Noise.mu = Noise.alpha*Noise.beta;
    Noise.cov = Noise.alpha*(Noise.beta^2);

    NoiseDS.type = 'NoiseDS';
    NoiseDS.ns_type = Noise.type;
    NoiseDS.tag = Noise.tag;
    NoiseDS.dim = Noise.dim;
    NoiseDS.alpha = Noise.alpha;
    NoiseDS.beta = Noise.beta;
    NoiseDS.mu = Noise.mu;
    NoiseDS.cov = Noise.cov;
    NoiseDS.sample = @sample_gamma;
    NoiseDS.likelihood = @likelihood_gamma;

%=============================================================================================
case 'gaussian'

    %-- check for and assign cov_type
    if isfield(ArgDS,'cov_type')
        if (ischar(ArgDS.cov_type) & stringmatch(ArgDS.cov_type,{'full','diag','sqrt','sqrt-diag'}))
            Noise.cov_type = ArgDS.cov_type;
        else
            error(' [ gennoiseds::gaussian ] Noise source cov_type not recognized or not a string.');
        end
    else
        warning(' [ gennoiseds::gaussian ] Covariance type field .cov_type not assigned!. Assuming default value, ''full''');
        Noise.cov_type = 'full';             % default cov_type
    end

    %-- assign noise source mean vector
    if ~isfield(ArgDS,'mu')
        Noise.mu = zeros(Noise.dim,1);    % default value
    else
        Noise.mu = ArgDS.mu;
    end

    if ~isfield(ArgDS,'cov')
      warning(' [ gennoiseds::gaussian ] Covariance field .cov not assigned!. Assuming default unity value.');
    end

    %-- assign rest of noise source structure
    switch (Noise.cov_type)

    %.............................................................................................
    case 'full'

        %-- assign noise source covariance structure
        if ~isfield(ArgDS,'cov')
            Noise.cov = eye(Noise.dim);    % default value
        elseif (isnumeric(ArgDS.cov) & (size(ArgDS.cov) == [Noise.dim Noise.dim]))
            Noise.cov = ArgDS.cov;
        else
            error([' [ gennoiseds::gaussian::full ] Noise source covariance matrix of incorrect ' ...
                   'dimensions or type.']);
        end

    %.............................................................................................
    case 'diag'

        %-- assign noise source covariance structure
        if ~isfield(ArgDS,'cov')
            Noise.cov = eye(Noise.dim);           % default value
        elseif (isnumeric(ArgDS.cov) & (size(ArgDS.cov) == [Noise.dim Noise.dim]))
            % check if covariance only has entries on the diagonal
            if (ArgDS.cov==diag(diag(ArgDS.cov)))
                Noise.cov = ArgDS.cov;              % assign covariance matrix
            else
                error([' [ gennoiseds::gaussian::diag ] Diagonal Guassian noise source cannot ' ...
                       'have non-zero off diagonal values in the covariance matrix.']);
            end
        else
            error([' [ gennoiseds::gaussian::diag ] Noise source covariance matrix of incorrect ' ...
                  'dimensions or type.']);
        end

    %.............................................................................................
    case 'sqrt'

        %-- assign noise source covariance structure
        if ~isfield(ArgDS,'cov')
            Noise.cov = eye(Noise.dim);    % default value
        elseif (isnumeric(ArgDS.cov) & (size(ArgDS.cov) == [Noise.dim Noise.dim]))
            Noise.cov = ArgDS.cov;
        else
            error([' [ gennoiseds::gaussian::sqrt ] Noise source covariance matrix of incorrect ' ...
                   'dimensions or type.']);
        end


    %.............................................................................................
    case 'sqrt-diag'

        %-- assign noise source covariance structure
        if ~isfield(ArgDS,'cov')
            Noise.cov = eye(Noise.dim);            % default value
        elseif (isnumeric(ArgDS.cov) & (size(ArgDS.cov) == [Noise.dim Noise.dim]))
            % check if covariance only has entries on the diagonal
            if (ArgDS.cov == diag(diag(ArgDS.cov)))
                Noise.cov = ArgDS.cov;               % assign covariance matrix
            else
                error([' [ gennoiseds::gaussian::sqrt-diag ] Diagonal Guassian noise source ' ...
                       'cannot have non-zero off diagonal values in the covariance matrix.']);
            end
        else
            error([' [ gennoiseds::gaussian::sqrt-diag ] Noise source covariance matrix of ' ...
                   'incorrect dimensions or type.']);
        end

    %.............................................................................................
    otherwise
        error(' [ gennoiseds::gaussian ]unknown noise source cov_type.');
    end


    % Restructure NoiseDS data structure
    NoiseDS.type = 'NoiseDS';
    NoiseDS.ns_type = Noise.type;
    NoiseDS.cov_type = Noise.cov_type;
    NoiseDS.tag = Noise.tag;
    NoiseDS.dim = Noise.dim;
    NoiseDS.mu = Noise.mu;
    NoiseDS.cov = Noise.cov;
    NoiseDS.sample = @gaussample;
    NoiseDS.update = @update_gaussian;
    NoiseDS.likelihood = @gauseval;




%===================================================================================================
case 'combo-gaussian'

    if (~isfield(ArgDS,'noiseSources') | ~iscell(ArgDS.noiseSources))
        error([' [ gennoiseds::combo-gaussian ] Sub noise source field (ArgDS.noiseSources) is ' ...
               'missing or is not a cell array.']);
    end

    Noise.N = length(ArgDS.noiseSources);    % number of sub noise sources
    if (Noise.N < 2)
        error([' [ gennoiseds::combo-gaussian ] A combo-Gaussian noise sources needs at ' ...
               'least 2 sub noise sources.']);
    end

    noisetype = ArgDS.noiseSources{1}.ns_type;
    cov_type = ArgDS.noiseSources{1}.cov_type;

    % Check cov_type type for correctness
    if ~(stringmatch(noisetype,{'gaussian','combo-gaussian'}) & ...
         stringmatch(cov_type,{'full','diag','sqrt','sqrt-diag'}))
        error(['[ gennoiseds::combo-gaussian ] A combination Gaussian noise source can ' ...
               'only have Gaussian sub noise sources.']);
    end

    % check for consistentency of cov_type of sub noise sources
    for k=1:Noise.N
        subNoise = ArgDS.noiseSources{k};
        if ~stringmatch(subNoise.cov_type,cov_type)
            error('[ gennoiseds::combo-gaussian ] Sub noise sources does not have consistent cov_types.\n Previous cov_type: %s     Current cov_type: %s', cov_type,subNoise.cov_type);
        end
    end

    Noise.cov_type = cov_type;                    % assign cov_type

    Noise.mu   = zeros(Noise.dim,1);            % setup mean vector
    Noise.cov  = zeros(Noise.dim);              % setup covariance matrix

    Noise.idxArr = zeros(Noise.N,2);            % buffer for beginning and ending indeces of sub noise
                                                % source entries in global mean and covariance

    % Extract sub noise source detail and build combo noise source
    dim = 0;
    for j=1:Noise.N,
        subNoise = ArgDS.noiseSources{j};       % copy j'th sub noise source
        dim = dim + subNoise.dim;               % add sub noise dimension to global dimension
        ind1 = dim-subNoise.dim+1;              % calculate beginning index value
        ind2 = dim;                             % calculate ending index value
        Noise.idxArr(j,:) = [ind1 ind2];        % store indeces in index array
        Noise.mu(ind1:ind2,1) = subNoise.mu;    % copy sub noise source mean
        Noise.cov(ind1:ind2,ind1:ind2) = subNoise.cov;  % copy sub noise covariance
    end
    if (Noise.dim ~= dim),
        error([' [ gennoiseds::combo-gaussian ] Combined noise vector dimension does ' ...
               'not agree with aggregate dimension of sub noise sources.']);
    end

    % Restructure NoiseDS data structure
    NoiseDS.type = 'NoiseDS';
    NoiseDS.ns_type = Noise.type;
    NoiseDS.cov_type = Noise.cov_type;
    NoiseDS.tag = Noise.tag;
    NoiseDS.N = Noise.N;
    NoiseDS.idxArr = Noise.idxArr;
    NoiseDS.dim = Noise.dim;
    NoiseDS.mu = Noise.mu;
    NoiseDS.cov = Noise.cov;
    NoiseDS.noiseSources = ArgDS.noiseSources;
    NoiseDS.sample = @sample_combo_gaussian;
    NoiseDS.update = @update_combo_gaussian;
    NoiseDS.likelihood = @likelihood_combo_gaussian;

%==================================================================================================
case 'combo'

    if (~isfield(ArgDS,'noiseSources') | ~iscell(ArgDS.noiseSources))
        error([' [ gennoiseds::combo ] Sub noise source field (ArgDS.noiseSources) ' ...
               'is missing or is not a cell array.']);
    end

    Noise.N = length(ArgDS.noiseSources);     % number of sub noise sources
    if (Noise.N < 2)
      error('[ gennoiseds::combo-gaussian ] A combo-Gaussian noise sources needs at least 2 sub noise sources.');
    end

    Noise.noiseSources = ArgDS.noiseSources;  % copy cell-array of sub noise sources

    Gaussian_flag = 1;                        % flag to check for the possibility of a Gaussian combination

    Noise.idxArr = zeros(Noise.N,2);          % buffer for beginning and ending indeces of sub noise
                                              % source entries in global mean and covariance

    dim = 0;
    for j=1:Noise.N,
        subNoise = Noise.noiseSources{j};
        dim = dim + subNoise.dim;
        ind1 = dim-subNoise.dim+1;
        ind2 = dim;
        Noise.idxArr(j,:) = [ind1 ind2];
        if ~stringmatch(subNoise.ns_type,{'gaussian','combo-gaussian'})
            Gaussian_flag = 0;
        end
    end

    if (Noise.dim ~= dim),
        error([' [ gennoiseds::combo-gaussian ] Combined noise vector dimension does ' ...
               'not agree with aggregate dimension of sub noise sources.']);
    end

    % A Gaussian combination should always be constructed if  the underlying sub
    % noise sources are all Gaussian
    if Gaussian_flag
        Arg.type = 'combo-gaussian';
        Arg.tag  = Noise.tag;
        Arg.dim  = Noise.dim;
        Arg.noiseSources = Noise.noiseSources;
        NoiseDS = gennoiseds(Arg);
    else
        % Restructure NoiseDS data structure
        NoiseDS.type = 'NoiseDS';
        NoiseDS.ns_type = Noise.type;
        NoiseDS.cov_type = Noise.cov_type;
        NoiseDS.tag = Noise.tag;
        NoiseDS.N = Noise.N;
        NoiseDS.idxArr = Noise.idxArr;
        NoiseDS.dim = Noise.dim;
        NoiseDS.mu = Noise.mu;
        NoiseDS.cov = Noise.cov;
        NoiseDS.noiseSources = Noise.noiseSources;
        NoiseDS.sample = @sample_combo;
        NoiseDS.update = @update_combo;
        NoiseDS.likelihood = @likelihood_combo;
    end


%===============================================================================================
%--- Gaussian Mixture Model
case 'gmm'

    if isfield(ArgDS,'M')
        Noise.M = ArgDS.M;
    else
        error(' [ gennoiseds::gmm ] Number of mixture components not specified.');
    end

    %-- assign noise source mean vector
    if ~isfield(ArgDS,'mu')
        Noise.mu = zeros(Noise.dim,Noise.M);    % default value
    else
        if (size(ArgDS.mu)==[Noise.dim Noise.M])
            Noise.mu = ArgDS.mu;
        else
            error(' [ gennoiseds::gmm ] Centroid mean dimension error.');
        end
    end

    %-- check for and assign cov_type
    if isfield(ArgDS,'cov_type')
        Noise.cov_type = ArgDS.cov_type;
    else
        warning(' [ gennoiseds::gmm ] Covariance type field .cov_type not assigned!. Assuming default value, ''full''');
        Noise.cov_type = 'full';             % default cov_type
    end

    %-- assign mixing weights
    if ~isfield(ArgDS,'weights')
        Noise.weights = (1/Noise.M)*ones(1,Noise.M);
    else
        if (length(ArgDS.weights)==Noise.M)
            Noise.weights = ArgDS.weights/(sum(ArgDS.weights));   % assign normalized weights
        else
            error(' [ gennoiseds::gmm ] Incorrect number of mixing weights (priors).');
        end
    end


    %-- assign rest of noise source structure
    switch (Noise.cov_type)

    %.............................................................................................
    case {'full','diag'}

        %-- assign noise source covariance structure
        if ~isfield(ArgDS,'cov')
            Noise.cov = repmat(eye(Noise.dim),[1 1 Noise.M]);    % default value
            warning(' [ gennoiseds::gmm ] Covariance field .cov not assigned!. Assuming default unity value.');
        elseif ((Noise.M == 1) & (size(ArgDS.cov) == [Noise.dim Noise.dim])) | ...
               ((Noise.M  > 1) & (size(ArgDS.cov) == [Noise.dim Noise.dim Noise.M]))
            Noise.cov = ArgDS.cov;
        else
            error([' [ gennoiseds::gmm::full ] Noise source covariance matrix has incorrect dimensions.']);
        end

    %.............................................................................................
    case {'sqrt','sqrt-diag'}

        %-- assign noise source covariance structure
        if ~isfield(ArgDS,'cov')
            Noise.cov = repmat(eye(Noise.dim),[1 1 Noise.M]);    % default value
            warning(' [ gennoiseds::gmm ] Covariance field .cov not assigned!. Assuming default unity value.');
        elseif ((Noise.M == 1) & (size(ArgDS.cov) == [Noise.dim Noise.dim])) | ...
               ((Noise.M  > 1) & (size(ArgDS.cov) == [Noise.dim Noise.dim Noise.M]))
            Noise.cov = ArgDS.cov;
        else
            error([' [ gennoiseds::gmm::sqrt ] Noise source covariance matrix has incorrect dimensions.']);
        end

    %.............................................................................................
    otherwise
        error(' [ gennoiseds::gmm ]unknown noise source cov_type.');
    end


    % Restructure NoiseDS data structure
    NoiseDS.type = 'NoiseDS';
    NoiseDS.ns_type = Noise.type;
    NoiseDS.cov_type = Noise.cov_type;
    NoiseDS.tag = Noise.tag;
    NoiseDS.dim = Noise.dim;
    NoiseDS.M = Noise.M;
    NoiseDS.weights = Noise.weights;
    NoiseDS.mu = Noise.mu;
    NoiseDS.cov = Noise.cov;
    NoiseDS.sample = @gmmsample;
    NoiseDS.likelihood = @likelihood_gmm;




%===============================================================================================

otherwise
  error([' [ gennoiseds ] Noise type ' type ' not supported.']);
end

NoiseDS.adaptMethod = [];



%***********************************************************************************************
%***                                                                                         ***
%***                               SUB FUNCTION BLOCK                                        ***
%***                                                                                         ***
%***********************************************************************************************

%===============================================================================================
function NoiseDS = update_gaussian(NoiseDS)

% Updates a Gaussian noise source
%
% This function is only a placeholder here since Gaussian noise sources are completely updated, i.e.
% mean and covariance are set externally, i.e. there are no INTERNAL structure to update.


%===============================================================================================
function noise = sample_gaussian(NoiseDS, N)

% Generate N samples of a noise source specified by the NoiseDS data structure

    switch NoiseDS.cov_type

    case 'full'
        A = chol(NoiseDS.cov)';
    case 'diag'
        A = diag(sqrt(diag(NoiseDS.cov)));
    case 'sqrt'
        A = NoiseDS.cov;
    case 'sqrt-diag'
        A = NoiseDS.cov;
    otherwise
        error(' [ sample_gaussian ] Unknown cov_type.');
    end

    noise = A*randn(NoiseDS.dim,N) + cvecrep(NoiseDS.mu,N);


%===============================================================================================
function llh = likelihood_combo_gaussian(NoiseDS, noise, idxVec)

% Calculates the likelihood of sample 'noise', given the noise model NoiseDS. If the optional index
% vector 'idxVec' is specified, only those sub-noise sources are used. The 'noise' vector's dimension should
% concur with the implied total dimensionality of 'idxVec'

    if (nargin == 2)
      idxVec = 1:NoiseDS.N;
    end

    numNS = length(idxVec);

    [dim,nov] = size(noise);
    idxArr = NoiseDS.idxArr;        % copy beginning/ending index array
    llh = ones(1,nov);

    % ... for each noise source
    for j=1:numNS,

        idx1 = idxArr(idxVec(j),1);
        idx2 = idxArr(idxVec(j),2);

        idxRange = idx1:idx2;

        dim = idx2-idx1+1;

        switch NoiseDS.cov_type

        case 'full'
            D  = det(NoiseDS.cov(idxRange,idxRange));
            iP = inv(NoiseDS.cov(idxRange,idxRange));
        case 'diag'
            D = prod(diag(NoiseDS.cov(idxRange,idxRange)));
            iP = diag(1./diag(NoiseDS.cov(idxRange,idxRange)));
        case 'sqrt'
            D = det(NoiseDS.cov(idxRange,idxRange))^2;
            iS = inv(NoiseDS.cov(idxRange,idxRange));
            iP = iS'*iS;
        case 'sqrt-diag'
            D = prod(diag(NoiseDS.cov(idxRange,idxRange)))^2;
            iP = diag(1./(diag(NoiseDS.cov(idxRange,idxRange)).^2));
        otherwise
            error(' [ likelihood_gaussian ]unknown cov_type.');
        end

        X = noise - cvecrep(NoiseDS.mu(idxRange),nov);
        q = 1/sqrt((2*pi)^dim * D);

        llh = llh .* (q * exp(-0.5*diag(X'*iP*X)'));

    end

    llh = llh + 1e-99; % needed to avoid 0 likelihood (cause ill conditioning)


%===============================================================================================
function llh = likelihood_combo(NoiseDS, noise, idxVec)

% Calculates the likelihood of sample 'noise', given the noise model NoiseDS.
% 'idxVec' is an optional index vector that can be used to indicate which of the N sub-noise sources should be used.
% to calculate the likelihood... this also requires 'noise' to have the same dimension of the relevant sub-noise source.

    if (nargin == 2)
      idxVec = 1:NoiseDS.N;
    end

    numNS = length(idxVec);

    [dim,nov] = size(noise);
    idxArr = NoiseDS.idxArr;        % copy beginning/ending index array
    llh = ones(1,nov);

    % ... for each noise source
    for j=1:numNS,

        idx1 = idxArr(idxVec(j),1);
        idx2 = idxArr(idxVec(j),2);
        subNoiseDS = NoiseDS.noiseSources{idxVec(j)};
        llh = llh .* subNoiseDS.likelihood( subNoiseDS, noise(idx1:idx2,:));

    end

    llh = llh + 1e-99; % needed to avoid 0 likelihood (cause ill conditioning)


%===============================================================================================
function noise = sample_combo(NoiseDS, N)

%  Generate N samples of a noise source specified by the NoiseDS data structure

    noise=zeros(NoiseDS.dim,N);     % setup noise sample output buffer

    idxArr = NoiseDS.idxArr;        % copy beginning/ending index array

    % ... for each noise source
    for j=1:NoiseDS.N
        subNoiseDS = NoiseDS.noiseSources{j};
        noise(idxArr(j,1):idxArr(j,2),:) = subNoiseDS.sample( subNoiseDS, N);
    end



%===============================================================================================
function noise = sample_combo_gaussian(NoiseDS, N)

%  Generate N samples of a noise source specified by the NoiseDS data structure

    noise=cvecrep(NoiseDS.mu,N);    % setup noise sample output buffer

    idxArr = NoiseDS.idxArr;        % copy beginning/ending index array

    num = NoiseDS.N;

    for j=1:num
        ind1 = idxArr(j,1);
        ind2 = idxArr(j,2);
        switch NoiseDS.cov_type
        case 'full'
            A = chol(NoiseDS.cov(ind1:ind2,ind1:ind2))';
        case 'diag'
            A = diag(sqrt(diag(NoiseDS.cov(ind1:ind2,ind1:ind2))));
        case 'sqrt'
            A = NoiseDS.cov(ind1:ind2,ind1:ind2);
        case 'sqrt-diag'
            A = NoiseDS.cov(ind1:ind2,ind1:ind2);
        otherwise
            error(' [ sample_gaussian ]unknown cov_type.');
        end
        noise(ind1:ind2,:) = noise(ind1:ind2,:) + A*randn(ind2-ind1+1,N);
    end



%===============================================================================================
function NoiseDS = update_combo_gaussian(NoiseDS)

% Updates a 'combination Gaussian' noise source which has N Gaussian sub noise sources. The global mean and covariance
% is updated externally and then this function is called to update the internal sub-noise source structure.
%

idxArr = NoiseDS.idxArr;

for j=1:NoiseDS.N,
    ind1 = idxArr(j,1);
    ind2 = idxArr(j,2);
    idxRange = ind1:ind2;
    NoiseDS.noiseSources{j}.mu = NoiseDS.mu(idxRange,1);
    NoiseDS.noiseSources{j}.cov = NoiseDS.cov(idxRange,idxRange);
end

%===============================================================================================
function noise = sample_gamma(NoiseDS, N)

% Generate N samples of a noise source specified by the NoiseDS data structure

    alpha = NoiseDS.alpha;
    beta = NoiseDS.beta;

    if (alpha==1)
        noise = -log(1-rand(1,N))*beta;
        return
    end

    flag=0;

    if (alpha<1)
        flag=1;
        alpha=alpha+1;
    end

    gamma=alpha-1;
    eta=sqrt(2.0*alpha-1.0);
    c=.5-atan(gamma/eta)/pi;

    y(N)=0;

    for k=1:N,
        aux=-.5;
        while(aux<0)
            y(k)=-.5;
            while(y(k)<=0)
                u=rand(1,1);
                y(k) = gamma + eta * tan(pi*(u-c)+c-.5);
            end
            v=-log(rand(1,1));
            aux=v+log(1.0+((y(k)-gamma)/eta)^2)+gamma*log(y(k)/gamma)-y(k)+gamma;
        end
    end

    if (flag==1)
        noise = y.*beta.*(rand(1,N)).^(1.0/(alpha-1));
    else
        noise = y.*beta;
    end



%===============================================================================================
function llh = likelihood_gamma(NoiseDS, noise)

% Calculates the likelihood of sample 'noise', given the noise model NoiseDS.

    llh = noise.^(NoiseDS.alpha-1) .* exp((-1/NoiseDS.beta)*noise);

%===============================================================================================
function llh = likelihood_gmm(NoiseDS, noise)

% Calculates the likelihood of sample 'noise', given the noise model NoiseDS.

[prior,likelihood,evidence] = gmmprobability(NoiseDS, noise);

llh = evidence;


