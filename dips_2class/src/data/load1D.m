function [data] = load1D(dirname)

%------------------------------------------------------------------------%
% Load time-series data from ABIDE I Preprocessed 
% These subjects have to pass a quality control step ahead.
%------------------------------------------------------------------------%

% remove empty file or time-series
system(['find ',dirname,' -type f -empty -delete']);

% import time-series from text files
fstruct = dir([dirname, '/*.1D']);
fnum    = length(fstruct);
data    = {};

for i = 1:fnum

	% load time-series of each subject
    disp(['reading: time-series ', fstruct(i).name]);
    d = importdata([dirname, '/', fstruct(i).name]);
    subj.data   = d.data;
    subj.name   = fstruct(i).name;

    % split file names to extract subject information
    s = strsplit(fstruct(i).name, {'_','.'});
    subj.site   = s{6};    	% scanning sties
    subj.id     = s{end-3};    % subject IDs
    subj.roi    = s{9};   		% atlas, parcellation
    
    data{end+1} = subj;

end