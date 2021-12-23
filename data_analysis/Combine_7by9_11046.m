% Combine_7by9__11046
%
% High level script to combine 7x9 data across sessions
%
% This version for the 7x9 pixel spots, for
% subject 11046.

%% Clear
clear; close all;

%% Parameters
%
% Where's the analyzed data and precomputed computational observer output?
theSubject = '11046';
psychoProject = 'AOPsychophysics';
psychoBaseDir = getpref(psychoProject,'analysisDir');
dataFilename = '11046_incDecFits_Aggregated.mat';

% Computational observer parameters
compProject = 'AOCompObserver';
compBaseDir = getpref(compProject,'analysisDir');
computationalName = '7_9_0';
defocusDiopters = 0.10;
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
titleStr = '11046, 7x9';
theLim = 2;

% Specify sessions.
sessionNames = {'IncrDecr1', 'IncrDecr2', 'IncrDecr4', ...
    'IncrDecr5', 'IncrDecr5', 'IncrDecr5', 'IncrDecr5', 'IncrDecr5', ...
    'IncrDecr6', 'IncrDecr6', 'IncrDecr6', 'IncrDecr6', 'IncrDecr6', ...
    'IncrDecr6', 'IncrDecr6', 'IncrDecr6', 'IncrDecr6', 'IncrDecr6', ...
    'IncrDecr6', 'IncrDecr6', 'IncrDecr6', 'IncrDecr6', 'IncrDecr6', ...
    'IncrDecr6', 'IncrDecr6'};

idx = 1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20200131','Separation_1','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20210914','Size_1','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211026','Size_1','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211123','Size1_Sep0','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211123','Size1_Sep2','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211123','Size1_Sep4','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211123','Size1_Sep6','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211123','Size1_Sep8','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211217','Dir0_Sep0','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211217','Dir45_Sep0','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211217','Dir45_Sep4','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211217','Dir45_Sep8','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211217','Dir45_Sep12','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211217','Dir45_Sep16','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211217','Dir225_Sep0','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211217','Dir225_Sep4','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211217','Dir225_Sep8','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211217','Dir225_Sep12','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211217','Dir225_Sep16','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211217','Dir270_Sep0','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211217','Dir315_Sep0','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211217','Dir315_Sep4','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211217','Dir315_Sep8','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211217','Dir315_Sep12','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211217','Dir315_Sep16','notnorm_corrguess_norefl',dataFilename); idx = idx+1;

%% Call the driver
CombineEllDriver;
