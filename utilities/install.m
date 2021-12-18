mainpath = fileparts(mfilename('fullpath'));
addpath(mainpath)
subfolders = {'initialization', 'kinematics', 'tests', 'correction', 'geodetic', ...
    'plotters', 'datastructure', 'imageproc', 'io', 'propagation'};
for i=1:length(subfolders)
    subpath = fullfile(mainpath, subfolders{i});
    fprintf('Add %s', subpath);
    addpath(subpath);
end
