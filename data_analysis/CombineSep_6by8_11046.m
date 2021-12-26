% CombineSep_6by8__11046
%
% High level script to combine 6x8 data across sessions
%
% This version for the 6x8 pixel spots, for
% subject 11046.
%
% This is still using 7x9 data, just set it up to allow looking
% at 6x8 computations.

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
computationalName = '6_8';
defocusDiopters = 0.05;
pupilDiam = 7;
compFilename = sprintf('%s_%s_D%s_P%d_ContourAnalysis.mat', ...
    theSubject,computationalName,num2str(round(1000*defocusDiopters)),pupilDiam);
compSeparations = [0 2 3 4 6 8 12 16];

% Normalization, etc.
SESSION_NORMALIZE = true;
INCRDECR_NORMALIZE = true;
ANGLE_AVERAGE = true;
REFLECT = false;
FIT_FIRSTTHIRDONLY = true;
angleTolerance = 1;
PFInputDir = 'notnorm_notcorrguess_norefl';

% Degrees per pixel
degsPerPixel = 1/415;
minPerPixel = 60*degsPerPixel;

% Other params
theLim = 3;

% Data directory infor and associated condition information filled in by
% hand
%
% Session numbers is an indicator variable so we can normalize within
% session if desired.
sessionNames = {'IncrDecr1', 'IncrDecr2', 'IncrDecr4', ...
    'IncrDecr5', 'IncrDecr5', 'IncrDecr5', 'IncrDecr5', 'IncrDecr5', ...
    'IncrDecr6', 'IncrDecr6', 'IncrDecr6', 'IncrDecr6', 'IncrDecr6', ...
    'IncrDecr6', 'IncrDecr6', 'IncrDecr6', 'IncrDecr6', 'IncrDecr6', ...
    'IncrDecr6', 'IncrDecr6', 'IncrDecr6', 'IncrDecr6', 'IncrDecr6', ...
    'IncrDecr6', 'IncrDecr6'};
runNames = {'Separation_1', 'Size_1', 'Size_1', 'Size1_Sep0', 'Size1_Sep2', 'Size1_Sep4', ...
    'Size1_Sep6', 'Size1_Sep8', 'Dir0_Sep0', 'Dir45_Sep0', 'Dir45_Sep4', 'Dir45_Sep8', ...
    'Dir45_Sep12', 'Dir45_Sep16', 'Dir225_Sep0', 'Dir225_Sep4', 'Dir225_Sep8', 'Dir225_Sep12', ...
    'Dir225_Sep16', 'Dir270_Sep0', 'Dir315_Sep0', 'Dir315_Sep4', 'Dir315_Sep8', 'Dir315_Sep12', ...
    'Dir315_Sep16'};
dateNames = {'20200131', '20210914', '20211026', '20211123', '20211123', '20211123', ...
    '20211123', '20211123', '20211217', '20211217', '20211217', '20211217', '20211217', ...  
    '20211217', '20211217', '20211217', '20211217', '20211217', '20211217', '20211217', ...   
    '20211217', '20211217', '20211217', '20211217', '20211217'};
sessionNumbers = [1 2 3 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5];
stimSeparationPixelsCheck = [0 0 0 0 2 4 6 8 0 0 4 8 12 16 0 4 8 12 16 0 0 4 8 12 16];
stimHeightCheck = 7*ones(size(stimSeparationPixelsCheck));
stimWidthCheck = 9*ones(size(stimSeparationPixelsCheck));

% Hand specification of which stimlus angles to plot separation series for
anglesToAnalyze = [45 225 315];

% Collect up the data file names
for idx = 1:length(sessionNames)
    theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,dateNames{idx},runNames{idx},PFInputDir,dataFilename);
end

%% Call the driver
CombineSepDriver;
