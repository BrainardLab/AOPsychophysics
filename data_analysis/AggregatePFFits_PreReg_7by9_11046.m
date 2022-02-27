% AggregatePFFits_PreReg_7by9_11046
%
% The initial PF fits constrain PF slope within condition, but sometimes
% there isn't enough data to do that, and it makes sense to leverage more
% data and get the overall fit with constrained PF slope. This script does
% the aggregation.
%
% This version for 7by9 stimuli for 11046, for pre-registered data.
%
% RunAllInitialPFFits must have been run first, to put the data in the form
% where this program can use it.

% History:
%    12/04/21  dhb  Wrote it.
%    12/20/21  dhb  Add some data checks, save more data.

%% Clear
clear; close all;

%% Parameters
%
% Where's are the initial PF fits?
psychoProject = 'AOPsychophysics';
psychoBaseDir = getpref(psychoProject,'analysisDir');

%% 11046, 7x9 data
theSubject = '11046';
theIndicator = 'incDecFits';

%% PF initial guess for grid search
initialParams = [];
initialAlphas = [];
initialBetas = [5 7 10 12 15];

%% Specify PF input
%
% Also set correction for guessing.
CORR_GUESSING = true;
PFInputDir = 'notnorm_notcorrguess_norefl';

%% Criterion percent seen for threshold
thresholdCriterion = 0.7;

% Data directory infor and associated condition information filled in by
% hand
%
% Session numbers is an indicator variable so we can normalize within
% session if desired.
sessionNames = {'IncrDecr8_7x9', ...
    'IncrDecr8_7x9', 'IncrDecr8_7x9', 'IncrDecr8_7x9', 'IncrDecr8_7x9', 'IncrDecr8_7x9', 'IncrDecr8_7x9', 'IncrDecr8_7x9', ...
    'IncrDecr8_7x9', 'IncrDecr8_7x9', 'IncrDecr8_7x9', 'IncrDecr8_7x9', 'IncrDecr8_7x9', 'IncrDecr8_7x9', 'IncrDecr8_7x9', ...
    'IncrDecr8_7x9', ...
    'IncrDecr8_7x9', 'IncrDecr8_7x9', 'IncrDecr8_7x9', 'IncrDecr8_7x9', 'IncrDecr8_7x9', 'IncrDecr8_7x9', 'IncrDecr8_7x9', ...
    ...
    'IncrDecr8_7x9', ...
    'IncrDecr8_7x9', 'IncrDecr8_7x9', 'IncrDecr8_7x9', 'IncrDecr8_7x9', 'IncrDecr8_7x9', 'IncrDecr8_7x9', 'IncrDecr8_7x9', ...
    'IncrDecr8_7x9', 'IncrDecr8_7x9', 'IncrDecr8_7x9', 'IncrDecr8_7x9', 'IncrDecr8_7x9', 'IncrDecr8_7x9', 'IncrDecr8_7x9', ...
    'IncrDecr8_7x9', ...
    'IncrDecr8_7x9', 'IncrDecr8_7x9', 'IncrDecr8_7x9', 'IncrDecr8_7x9', 'IncrDecr8_7x9', 'IncrDecr8_7x9', 'IncrDecr8_7x9'};
runNames = {'Dir0_Sep0', ...
    'Dir45_Sep0', 'Dir45_Sep2', 'Dir45_Sep4', 'Dir45_Sep7', 'Dir45_Sep13', 'Dir45_Sep24', 'Dir45_Sep45', ...
    'Dir225_Sep0', 'Dir225_Sep2', 'Dir225_Sep4', 'Dir225_Sep7', 'Dir225_Sep13', 'Dir225_Sep24', 'Dir225_Sep45', ...
    'Dir270_Sep0', ...
    'Dir315_Sep0', 'Dir315_Sep2', 'Dir315_Sep4', 'Dir315_Sep7', 'Dir315_Sep13', 'Dir315_Sep24', 'Dir315_Sep45', ...
    ...
    'Dir0_Sep0', ...
    'Dir45_Sep0', 'Dir45_Sep2', 'Dir45_Sep4', 'Dir45_Sep7', 'Dir45_Sep13', 'Dir45_Sep24', 'Dir45_Sep45', ...
    'Dir225_Sep0', 'Dir225_Sep2', 'Dir225_Sep4', 'Dir225_Sep7', 'Dir225_Sep13', 'Dir225_Sep24', 'Dir225_Sep45', ...
    'Dir270_Sep0', ...
    'Dir315_Sep0', 'Dir315_Sep2', 'Dir315_Sep4', 'Dir315_Sep7', 'Dir315_Sep13', 'Dir315_Sep24', 'Dir315_Sep45'};
dateNames = {'20220211', ...
    '20220211', '20220211', '20220211', '20220211', '20220211', '20220211', '20220211', ...
    '20220211', '20220211', '20220211', '20220211', '20220211', '20220211', '20220211', ...
    '20220211', ...
    '20220211', '20220211', '20220211', '20220211', '20220211', '20220211', '20220211', ...
    ...
    '20220223', ...
    '20220223', '20220223', '20220223', '20220223', '20220223', '20220223', '20220223', ...
    '20220223', '20220223', '20220223', '20220223', '20220223', '20220223', '20220223', ...
    '20220223', ...
    '20220223', '20220223', '20220223', '20220223', '20220223', '20220223', '20220223'};
sessionNumbers = [1*ones(size(sessionNames)) 2*ones(size(sessionNames))];
stimSeparationPixelsCheckRaw = [0 0 2 4 7 13 24 45 0 2 4 7 13 24 45 0 0 2 4 7 13 24 45];
stimSeparationPixelsCheck = [stimSeparationPixelsCheckRaw stimSeparationPixelsCheckRaw];
stimAngleCheckRaw = [0 45 45 45 45 45 45 45 225 225 225 225 225 225 225 270 315 315 315 315 315 315 315];
stimAngleCheck = [stimAngleCheckRaw stimAngleCheckRaw];
stimHeightCheck = 7*ones(size(stimSeparationPixelsCheck));
stimWidthCheck = 9*ones(size(stimSeparationPixelsCheck));

% Collect up the data file names
for idx = 1:length(sessionNames)
    theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,dateNames{idx},runNames{idx},PFInputDir,sprintf('%s_%s_ConstrainedSlope.mat',theSubject,theIndicator));
end

%% Run the driver
AggregatePFFitsDriver;