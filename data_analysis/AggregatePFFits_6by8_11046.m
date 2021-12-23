% AggregatePFFits_6by8_11046
%
% The initial PF fits constrain PF slope within condition, but sometimes
% there isn't enough data to do that, and it makes sense to leverage more
% data and get the overall fit with constrained PF slope. This script does
% the aggregation.
%
% This version for 6by8 stimuli for 11046. That is the current view of the
% right set of data to lock the PF slope for.
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

%% 11046, 6x8 data
theSubject = '11046';
theIndicator = 'incDecFits';

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
sessionNames = {'IncrDecr4', 'IncrDecr5', 'IncrDecr5'};
runNames = {'Size_4', 'Size2_Sep0', 'Size2_Sep4'};
dateNames = {'20211026', '20211123', '20211123'};
sessionNumbers = [1 2 2];
stimSeparationPixelsCheck = [0 0 4];
stimHeightCheck = 6*ones(size(stimSeparationPixelsCheck));
stimWidthCheck = 8*ones(size(stimSeparationPixelsCheck));

% PF initial guess for grid search
initialParams = [];
initialAlphas = [0.2 0 -0.2];
initialBetas = [5 10 15];

% Collect up the data file names
for idx = 1:length(sessionNames)
    theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,dateNames{idx},runNames{idx},PFInputDir,sprintf('%s_%s_ConstrainedSlope.mat',theSubject,theIndicator));
end

%% Run the driver
AggregatePFFitsDriver;