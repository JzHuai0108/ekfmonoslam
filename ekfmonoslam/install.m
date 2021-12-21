mainpath = fileparts(mfilename('fullpath'));
addpath(mainpath)
subfolders = {'initialization', 'kinematics', 'tests', 'correction', 'ekf', 'geodetic', ...
    'plotters', 'datastructure', 'imageproc', 'io', 'propagation'};
for i=1:length(subfolders)
    subpath = fullfile(mainpath, subfolders{i});
    fprintf('Add %s\n', subpath);
    addpath(subpath);
end
