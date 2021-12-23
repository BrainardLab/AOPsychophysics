% AggregatePFFits_7by9_11046
%
% The initial PF fits constrain PF slope within condition, but sometimes
% there isn't enough data to do that, and it makes sense to leverage more
% data and get the overall fit with constrained PF slope. This script does
% the aggregation.
%
% This version for 7by9 stimuli for 11046. That is the current view of the
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

%% 11046, 7x9 data
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

% Collect up the data file names
for idx = 1:length(sessionNames)
    theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,dateNames{idx},runNames{idx},PFInputDir,sprintf('%s_%s_ConstrainedSlope.mat',theSubject,theIndicator));
end

%% Run the driver
AggregatePFFitsDriver;