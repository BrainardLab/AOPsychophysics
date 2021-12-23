% Combine_7by9_11002
%
% High level script to combine 7x9 data across sessions
%
% This version for the 7x9 pixel spots, for
% subject 11002.

%% Clear
clear; close all;

%% Parameters
%
% Where's the analyzed data and precomputed computational observer output?
theSubject = '11002';
psychoProject = 'AOPsychophysics';
psychoBaseDir = getpref(psychoProject,'analysisDir');
dataFilename = '11002_incDecFits_Aggregated.mat';

% Computational observer parameters
compProject = 'AOCompObserver';
compBaseDir = getpref(compProject,'analysisDir');
computationalName = '7_9_0';
defocusDiopters = 0.05;
pupilDiam = 7;
compFilename = sprintf('%s_%s_D%s_P%d_ContourAnalysis.mat', ...
    theSubject,computationalName,num2str(round(1000*defocusDiopters)),pupilDiam);

% Normalization, etc.
SESSION_NORMALIZE = true;
INCRDECR_NORMALIZE = true;
ANGLE_AVERAGE = true;
REFLECT = false;
angleTolerance = 1;

% Other params
titleStr = '11002, 7x9';
theLim = 2;

% Detailed data info
sessionNames = {'IncrDecr7_7x9', 'IncrDecr7_7x9', 'IncrDecr7_7x9', ...
    'IncrDecr7_7x9', 'IncrDecr7_7x9', 'IncrDecr7_7x9', 'IncrDecr7_7x9'};
runNames = {'Dir0_Sep0', 'Dir45_Sep0', 'Dir225_Sep0', ...
    'Dir270_Sep0', 'Dir293_Sep0', 'Dir315_Sep0', 'Dir338_Sep0'};
dateNames = {'20211220', '20211220', '20211220', ...
    '20211220', '20211220', '20211220', '20211220'};
sessionNumbers = [1 1 1 1 1 1 1];
stimSeparationPixelsCheck = [0 0 0 0 0 0 0];
stimHeightCheck = 7*ones(size(stimSeparationPixelsCheck));
stimWidthCheck = 9*ones(size(stimSeparationPixelsCheck));

% Collect up the data file names
for idx = 1:length(sessionNames)
    theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,dateNames{idx},runNames{idx},'notnorm_corrguess_norefl',dataFilename);
end

%% Call the driver
CombineEllDriver;
