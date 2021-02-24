% Script to combine data from inc/dec experiments and fit PFs to various
% axes in the 2D stimulus space

%% Housekeeping
clear; close all
baseProject = 'AOPsychophysics';
subProject = 'IncrDecr1';
dataBaseDir = getpref(baseProject,'dataDir');
analysisBaseDir = getpref(baseProject,'analysisDir');

%% Load in the data files (update directories to wherever these data live on your machine)
%
% Select subject
%   subj = '11043'; % WST
%   subj = '11046'; % DHB
subj = '11046';
switch (subj)
    case {'11043' '11046'}
        dataDate = '20200131';
        fprintf('Subject ID: %s\n', subj);
    otherwise
        error('Specified subject number invalid')
end

%% Normalize?
%
% To normalize, or not to normalize? Select false to work with the raw
% modulations, true to normalize.
normFlag = false;

% Set up directories
dataDir = fullfile(dataBaseDir,subProject,subj,dataDate,'Separation_1');
analysisDir = fullfile(analysisBaseDir,subProject,subj,dataDate,'Separation_1');
if (~exist(analysisDir,'dir'))
    mkdir(analysisDir);
end


% Get list of .mat files for the data we're analyzing
dataFiles = dir(fullfile(dataDir,'*DetectionData.mat'));
if isempty(dataFiles)
    error('No data files found');
end

%% Get the data
stimulusModulations = [];
responseVector = [];
for fileNum = 1:length(dataFiles)
    % Load individual data files
    expDataTemp = load(fullfile(dataFiles(fileNum).folder, dataFiles(fileNum).name));
    
    % Add to output matrix/vector
    stimulusModulations = [stimulusModulations; expDataTemp.testSeq(:,1:2)-expDataTemp.CFG.backgroundIntensity];
    responseVector = [responseVector; expDataTemp.YesNoResponseMatrix]; %#ok<*AGROW>
end

%% Normalize if desired
stimulusModulationsNormalized = zeros(size(stimulusModulations)); % Pre-allocate
if normFlag == 1    
    decNormFactor = max(abs(unique(stimulusModulations(stimulusModulations<0))));
    incNormFactor = max(unique(stimulusModulations(stimulusModulations>0)));
else
    decNormFactor = 1;
    incNormFactor = 1;
end
stimulusModulationsNormalized(stimulusModulations>0) = stimulusModulations(stimulusModulations>0)./incNormFactor;
stimulusModulationsNormalized(stimulusModulations<0) = stimulusModulations(stimulusModulations<0)./decNormFactor;

% Compute stimulus angle in 2-D space.  Just brute force wrt the
% experimental design.
stimAngles = zeros(length(stimulusModulations),1);
for n = 1:length(stimAngles)
    if stimulusModulationsNormalized(n,1) < 0
        stim1 = 'Dec';
    elseif stimulusModulationsNormalized(n,1) > 0
        stim1 = 'Inc';
    else
        stim1 = 'Blank';
    end
    
    if stimulusModulationsNormalized(n,2) < 0
        stim2 = 'Dec';
    elseif stimulusModulationsNormalized(n,2) > 0
        stim2 = 'Inc';
    else
        stim2 = 'Blank';
    end
    stimAngles(n) = round(atand(stimulusModulationsNormalized(n,2)./stimulusModulationsNormalized(n,1)));
    
    if stimAngles(n) == 45
        if strcmp(stim1,'Dec') && strcmp(stim2, 'Dec')
            stimAngles(n) = stimAngles(n)+180;
        end
    end
    
    if stimAngles(n) == -45
        stimAngles(n) = 360-45;
    end
    
    if strcmp(stim1, 'Dec') && strcmp(stim2, 'Blank')
        stimAngles(n) = 270;
    end
end

%% Compute the false positive rate
% Values of 3 and 1 are from the gamepad controller used in the experiment
yesVal = 3; % From gamepad
noVal = 1;
numCatchTrials = length(responseVector(mean(stimulusModulationsNormalized,2)==0));
catchResponses = responseVector(mean(stimulusModulationsNormalized,2)==0);
numCatchPos = length(catchResponses(catchResponses==yesVal));
falsePosProp = numCatchPos/numCatchTrials;
fprintf('False positive rate: %d out of %d; (%.2f percent)\n', numCatchPos, numCatchTrials, 100*(falsePosProp));

%% Compute PF for each modulation direction individually
%
% Convert absolute modulations into vector distance from the origin; can be
% thought of ans an absolute modulation "magnitude"
fprintf('Fitting modulation directions individually... ');
stimAngleList = unique(stimAngles(isfinite(stimAngles)));
modulationVectorLengths = sqrt(stimulusModulationsNormalized(:,1).^2 + stimulusModulationsNormalized(:,2).^2);

% Set up the PF fitting (requires Palamedes toolbox)
PF = @PAL_Logistic;
searchGrid = [log10(mean(modulationVectorLengths)) 5 falsePosProp, 0.01];
paramsFree = [1 1 0 0]; %[thresh slope guess lapse]; 0 = fixed; 1 = free

% Pre-allocate
paramsFitted_Individual = nan(length(stimAngleList),4);

% For PF evaluation
xEval = linspace(0,1,1000);

% Plot into this figure
figure; hold on;
set(gcf, 'Color', 'w', 'Units', 'inches', 'Position', [1 1 12 4]);

for angleNum = 1:length(stimAngleList)
    fitInd = find(stimAngles==stimAngleList(angleNum));
    
    angleResponses = responseVector(fitInd);
    angleModulations = modulationVectorLengths(fitInd);
    
    % Group data (not necessary but makes for easier plotting)
    stimLevels(angleNum,:) = [0.001 unique(angleModulations)']; %#ok<SAGROW>
    outOfNum(angleNum,:) = zeros(1,size(stimLevels,2)); %#ok<SAGROW>
    numPos(angleNum,:) = zeros(1,size(stimLevels,2)); %#ok<SAGROW>
    
    for j = 1:length(stimLevels)
        if j == 1 % Catch trial
            outOfNum(angleNum,j) = numCatchTrials;
            numPos(angleNum,j) = numCatchPos;
        else
            tempRespVector = angleResponses(angleModulations==stimLevels(angleNum,j));
            outOfNum(angleNum, j) = length(tempRespVector);
            numPos(angleNum, j) = length(tempRespVector(tempRespVector==yesVal));
        end
    end
    
    paramsFitted_Individual(angleNum,:) = PAL_PFML_Fit([-3 log10(stimLevels(angleNum,:))], [numCatchPos numPos(angleNum,:)], [numCatchTrials outOfNum(angleNum,:)], searchGrid, paramsFree, PF);
    
    % Add to plot
    ax1 = subplot(1,3,1); hold on;
    h = plot(log10(xEval), PF(paramsFitted_Individual(angleNum,:), log10(xEval)), '-', 'LineWidth', 2);
    plot(log10(stimLevels(angleNum,:)), numPos(angleNum,:)./outOfNum(angleNum,:), 's', ...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor', h.Color, 'MarkerSize', 10, 'HandleVisibility', 'off');
end

% Finish up plot
legend(num2str(stimAngleList), 'Location', 'NorthWest')
xlim([-1.5 0])
ylim([0 1]);
axis square
xlabel('Log10 modulation intensity (au)', 'FontSize', 14);
ylabel('Prop seen', 'FontSize', 14);
title('Slopes unconstrained; guess rate fixed');
fprintf('Done.\n');

%% Fit everything together with slopes and guess/lapse rates to be equal across stimulus angles
fprintf('Now fitting all angles with the same slope and guess rate... ');
paramsFitted_Multi = PAL_PFML_FitMultiple(log10(stimLevels), numPos, outOfNum, searchGrid, PF, 'slopes', 'constrained', 'guessrates', ...
    'fixed', 'lapserates', 'constrained', 'lapselimits', [0 0.05]);
fprintf('Done.\n');

% Plot
for angleNum = 1:length(stimAngleList)
    ax2 = subplot(1,3,2); hold on;
    h2 = plot(log10(xEval), PF(paramsFitted_Multi(angleNum,:), log10(xEval)), '-', 'LineWidth', 2);
    subplot(1,3,2), plot(log10(stimLevels(angleNum,:)), numPos(angleNum,:)./outOfNum(angleNum,:), 's', 'MarkerFaceColor', h2.Color, ...
        'MarkerEdgeColor', 'none', 'MarkerSize', 10, 'HandleVisibility', 'off');
end
legend(num2str(stimAngleList), 'Location', 'NorthWest')
xlim([-1.5 0])
ylim([0 1]);
title('Slopes constrained; guess rate fixed');
axis square;
xlabel('Log10 modulation intensity (au)', 'FontSize', 14);
ylabel('Prop seen', 'FontSize', 14);

%% Evaluate fit a specific prop seen and plot on 2-D modulation space
ax3 = subplot(1,3,3); hold on;
plot([0 0], [-1 1], 'k-', 'LineWidth', 1.5);
plot( [-1 1],[0 0], 'k-', 'LineWidth', 1.5);

% This is how I evaluate the PF and convert back to x,y coordinates for ellipse fitting
% First, round to the nearest .1 above the guess rate (code won't be happy if you
% try to evaluate the PF for prop seen below the lower asymptote)
propEvalStart = ceil(10*falsePosProp)/10; 
propSeen_Fit = propEvalStart:.1:.9;
modLevels_PF = nan(size(paramsFitted_Multi,1), length(propSeen_Fit));
for n = 1:size(paramsFitted_Multi,1)
    modLevels_PF(n,:) = 10.^PF(paramsFitted_Multi(n,:), propSeen_Fit, 'inv');
end

% Then, convert to x and y coordinates
xPlot_Fit = cosd(stimAngleList).*modLevels_PF;
yPlot_Fit = sind(stimAngleList).*modLevels_PF;

% Next, plot the coordinates, color-coded to prop seen
hold on; 
for n = 1:length(propSeen_Fit)
    subplot(1,3,3), plot(xPlot_Fit(:,n), yPlot_Fit(:,n), 'ks', 'MarkerSize', 10, 'MarkerFaceColor', [0 propSeen_Fit(n) 0]);
end

xlabel('Stimulus 1 modulation (normalized)', 'FontSize', 14);
ylabel('Stimulus 2 modulation (normalized)', 'FontSize', 14);
axLim = round(max(stimulusModulationsNormalized(:)),1); % Round to nearest 0.1
xlim([-axLim axLim]);
ylim([-axLim axLim]);
axis square;
set(gca, 'XTickMode', 'auto', 'YTickMode', 'auto')
box on; grid on;
title('PFs evaluated from center panel');

% Make a color bar and tweak its properties
cmap = zeros(length(0:0.1:1), 3);
cmap(:,2) = 0:.1:1;
c = colorbar;
colormap(cmap);
c.Label.String = 'Prop seen from PF fit';
c.Label.Rotation = 270;
c.Label.VerticalAlignment = 'bottom';
c.Label.FontSize = 12;

set(ax1, 'Position', [0.1 0.110 0.3347.*.66 0.8150])
set(ax2, 'Position', [0.4 0.110 0.3347.*.66 0.8150])
set(ax3, 'Position', [0.7 0.110 0.3347*.66 0.8150])

%% Save fit data to mat file
save(fullfile(analysisDir,sprintf('%s_incDecFits_ConstrainedSlope.mat', subj)), 'stimAngleList', 'falsePosProp', 'paramsFitted_Individual', 'paramsFitted_Multi', 'PF');
print(gcf, fullfile(analysisDir,sprintf('%s_incDecFits_ConstrainedSlope.png', subj)), '-dpng2');
