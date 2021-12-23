% Combine_7by9_11043
%
% High level script to combine 7x9 data across sessions
%
% This version for the 7x9 pixel spots, for
% subject 11043.

%% Clear
clear; close all;

%% Parameters
%
% Where's the analyzed data and precomputed computational observer output?
theSubject = '11043';
psychoProject = 'AOPsychophysics';
psychoBaseDir = getpref(psychoProject,'analysisDir');
dataFilename = [theSubject '_incDecFits_Aggregated.mat'];

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
PFInputDir = 'notnorm_notcorrguess_norefl';

% Other params
titleStr = '11043, 7x9';
theLim = 3;

% Detailed data info
sessionNames = {'IncrDecr1'};
runNames = {'Separation_1'};
dateNames = {'20200131'};
sessionNumbers = [1];
stimSeparationPixelsCheck = [0];
stimHeightCheck = 7*ones(size(stimSeparationPixelsCheck));
stimWidthCheck = 9*ones(size(stimSeparationPixelsCheck));

% Collect up the data file names
for idx = 1:length(sessionNames)
    theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,dateNames{idx},runNames{idx},PFInputDir,dataFilename);
end

%% Call the driver
CombineEllDriver;
