mainpath = fileparts(mfilename('fullpath'));
addpath(mainpath)
subfolders = {'Calibration', 'Examples102', 'initialization', 'Common', 'INS', 'Smoother', ...
    'conscull', 'IMUModeling', 'OldVersions', 'PathGen', 'PDANav', 'TransFunctions'};
for i=1:length(subfolders)
    subpath = fullfile(mainpath, subfolders{i});
    fprintf('Add %s\n', subpath);
    addpath(subpath);
end
