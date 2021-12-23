% AggregatePFFits_7by9_11043
%
% The initial PF fits constrain PF slope within condition, but sometimes
% there isn't enough data to do that, and it makes sense to leverage more
% data and get the overall fit with constrained PF slope. This script does
% the aggregation.
%
% This version for 7by9 stimuli for 11043. That is the current view of the
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

%% 11043, 7x9 data
theSubject = '11043';
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
sessionNames = {'IncrDecr1'};
runNames = {'Separation_1'};
dateNames = {'20200131'};
sessionNumbers = [1];
stimSeparationPixelsCheck = [0];
stimHeightCheck = 7*ones(size(stimSeparationPixelsCheck));
stimWidthCheck = 9*ones(size(stimSeparationPixelsCheck));

% Collect up the data file names
for idx = 1:length(sessionNames)
    theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,dateNames{idx},runNames{idx},PFInputDir,sprintf('%s_%s_ConstrainedSlope.mat',theSubject,theIndicator));
end

%% Run the driver
AggregatePFFitsDriver;
