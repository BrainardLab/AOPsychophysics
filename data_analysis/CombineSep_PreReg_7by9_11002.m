% CombineSep_PreReg_7by9_11002
%
% High level script to combine 7x9 data across sessions
%
% This version for the 7x9 pixel spots, for
% subject 11002.

%% Clear
clear; close all;

%% Parameters
%
% Where's the analyzed data and precomputed computational observer output?
theSubject = '11002';
psychoProject = 'AOPsychophysics';
psychoBaseDir = getpref(psychoProject,'analysisDir');
dataFilename = [theSubject '_incDecFits_Aggregated.mat'];

% Computational observer parameters
compProject = 'AOCompObserver';
compBaseDir = getpref(compProject,'analysisDir');
computationalName = '7_9';
defocusDiopters = 0.05;
pupilDiam = 7;
compFilename = sprintf('%s_%s_D%s_P%d_ContourAnalysis.mat', ...
    theSubject,computationalName,num2str(round(1000*defocusDiopters)),pupilDiam);
compSeparations = [0 1 2 3 4 5 6 7 8 12 13 16 20 24 35 45];

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
theLim = 6;

% Specify conditions, etc.
DataSpec_PreReg_7by9_11002;

% Hand specification of which stimlus angles to plot separation series for
anglesToAnalyze = [45 225 315];

% Collect up the data file names
for idx = 1:length(sessionNames)
    theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,dateNames{idx},runNames{idx},PFInputDir,dataFilename);
end

%% Call the driver
CombineSepDriver;
