% AggregatePFFits_7by9_11002
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
theSubject = '11002';
theIndicator = 'incDecFits';

%% Criterion percent seen for threshold
thresholdCriterion = 0.7;

%% Data directory infor and associated condition information
% 
% Filled in by hand
%
% Session numbers is an indicator variable so we can normalize within
% session if desired.
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
    theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,dateNames{idx},runNames{idx},'notnorm_corrguess_norefl',sprintf('%s_%s_ConstrainedSlope.mat',theSubject,theIndicator));
end

%% Call into the driver
AggregatePFFitsDriver;