% CombineEll_5by7__11046
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
computationalName = '5_7_0';
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
sessionNames = {'IncrDecr3', 'IncrDecr4'};
runNames = {'Size_2', 'Size_2'};
dateNames = {'20211018', '20211026'};
sessionNumbers = [1 2];
stimSeparationPixelsCheck = [0 0];
stimHeightCheck = 5*ones(size(stimSeparationPixelsCheck));
stimWidthCheck = 7*ones(size(stimSeparationPixelsCheck));

% Collect up the data file names
for idx = 1:length(sessionNames)
    theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,dateNames{idx},runNames{idx},PFInputDir,dataFilename);
end

%% Call the driver
CombineEllDriver;
