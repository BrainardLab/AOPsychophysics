% CombineEll_6by8__11046
%
% High level script to combine 6x8 data across sessions
%
% This version for the 6x8 pixel spots, for
% subject 11046.

%% Clear
clear; close all;

%% Parameters
%
% Where's the analyzed data and precomputed computational observer output?
theSubject = '11046';
psychoProject = 'AOPsychophysics';
psychoBaseDir = getpref(psychoProject,'analysisDir');
dataFilename = [theSubject '_incDecFits_Aggregated.mat'];

% Computational observer parameters
compProject = 'AOCompObserver';
compBaseDir = getpref(compProject,'analysisDir');
computationalName = '6_8_0';
defocusDiopters = 0.05;
pupilDiam = 7;
compFilename = sprintf('%s_%s_D%s_P%d_ContourAnalysis.mat', ...
    theSubject,computationalName,num2str(round(1000*defocusDiopters)),pupilDiam);

% Normalization, etc.
SESSION_NORMALIZE = true;
INCRDECR_NORMALIZE = true;
ANGLE_AVERAGE = true;
REFLECT = false;
FIT_FIRSTTHIRDONLY = true;
angleTolerance = 1;
PFInputDir = 'notnorm_notcorrguess_norefl';

% Other params
theLim = 5;

% Data directory infor and associated condition information filled in by
% hand
%
% Session numbers is an indicator variable so we can normalize within
% session if desired.
sessionNames = {'IncrDecr4', 'IncrDecr5', 'IncrDecr5'};
runNames = {'Size_4', 'Size2_Sep0', 'Size2_Sep4'};
dateNames = {'20211026', '20211123', '20211123'};
sessionNumbers = [1 2 2];
stimSeparationPixelsCheck = [0 0 4];
stimHeightCheck = 6*ones(size(stimSeparationPixelsCheck));
stimWidthCheck = 8*ones(size(stimSeparationPixelsCheck));

% Collect up the data file names
for idx = 1:length(sessionNames)
    theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,dateNames{idx},runNames{idx},PFInputDir,dataFilename);
end

%% Call the driver
CombineEllDriver;
