function fitPfsToIncrDecrData(options)
% Combine data from inc/dec experiments and fit PFs
%
% Description:
%   Read in data from the AO inc/dec experiments and fit PFs.
%
%   Optional key/value pairs
%      'subj'       - String. Subject ID.  Default '11043';
%      'dataDate'   - String. Date data collected. Default '20200131'.
%      'subProject' - String. Subproject. Default 'IncrDecr1'.
%      'condition'  - String. Specify condition. Default 'Separation_1'.
%      'norm'       - Boolean. Normalize increment constrast by max used, and same for
%                     decrement contrast. Default false.
%      'corrGuess'  - Boolean. Correct for guessing using high threshold
%                     model. Default true.
%      'refl'       - Boolean. Treat stim1 and stim2 as symmetric, and reflect data
%                     Default false.
%
%   Data are fit with independent psychometric functions in each direction
%   as well as with fits constrained across directions.  
%
%   A little more work to make sure we're happy with how guess and lapse
%   rates are set/constrained in the fits would be good before this
%   analysis is finalized.

% History:
%   xx/xx/xx  wst  Wrote initial version.
%   03/05/21  dhb  Cleaning up and commenting after some earlier changes.

%% Pick up optional arguments 
arguments
    options.subj string = '11043';
    options.dataDate string = '20200131';
    options.subProject string = 'IncrDecr1';
    options.condition string = 'Separation_1';
    options.norm (1,1) logical = false;
    options.corrGuess (1,1) logical = true;
    options.refl (1,1) logical = false;
end

%% Housekeeping
close all
baseProject = 'AOPsychophysics';
subProject = options.subProject;

%% Load in the data files (update directories to wherever these data live on your machine)
switch (options.subj)
    case {'11043' '11046'}
        fprintf('Subject ID: %s\n', options.subj);
    otherwise
        error('Specified subject number invalid')
end

%% Normalize?
%
% To normalize, or not to normalize? Select false to work with the raw
% modulations, true to normalize.
if (options.norm)
    normStr = 'norm';
else
    normStr = 'notnorm';
end

%% Correct for guessing?
if (options.corrGuess)
    corrGuessStr = 'corrguess';
else
    corrGuessStr = 'notcorrguess';
end

%% Refelct data so that stim1 and stim2 are treated as symmetric?
if (options.refl)
    reflStr = 'refl';
else
    reflStr = 'norefl';
end

%% Set up directories
analysisSubDir = sprintf('%s_%s_%s',normStr,corrGuessStr,reflStr);
dataBaseDir = getpref(baseProject,'dataDir');
analysisBaseDir = getpref(baseProject,'analysisDir');
dataDir = fullfile(dataBaseDir,subProject,options.subj,options.dataDate,options.condition);
analysisDir = fullfile(analysisBaseDir,subProject,options.subj,options.dataDate,options.condition,analysisSubDir);
if (~exist(analysisDir,'dir'))
    mkdir(analysisDir);
end

%% Get list of .mat files for the data we're analyzing
dataFiles = dir(fullfile(dataDir,'*DetectionData.mat'));
if isempty(dataFiles)
    error('No data files found');
end

%% Get the data.
%
% First entry is top, second is bottom
stimulusContrasts = [];
responseVector = [];
for fileNum = 1:length(dataFiles)
    % Load individual data files
    expDataTemp = load(fullfile(dataFiles(fileNum).folder, dataFiles(fileNum).name));
    
    % Get stimulus contrast, without symmeterizing.  First entry of data is
    % top stimulus if the fifth column is 2, but the other way around if
    % fifth column is 4.
    stimulusContrastsIn = (expDataTemp.testSeq(:,1:2)-expDataTemp.CFG.backgroundIntensity)/expDataTemp.CFG.backgroundIntensity;
    flipIndex = find(expDataTemp.testSeq(:,5) == 4);
    stimulusContrastsIn(flipIndex,:) = stimulusContrastsIn(flipIndex,[2 1]);

    % Add to output matrix/vector
    stimulusContrasts = [stimulusContrasts; stimulusContrastsIn];
    responseVector = [responseVector; expDataTemp.YesNoResponseMatrix]; %#ok<*AGROW>
end

%% Normalize if desired
stimulusModulationsNormalized = zeros(size(stimulusContrasts)); % Pre-allocate
if options.norm == 1    
    decNormFactor = max(abs(unique(stimulusContrasts(stimulusContrasts<0))));
    incNormFactor = max(unique(stimulusContrasts(stimulusContrasts>0)));
else
    decNormFactor = 1;
    incNormFactor = 1;
end
stimulusModulationsNormalized(stimulusContrasts>0) = stimulusContrasts(stimulusContrasts>0)./incNormFactor;
stimulusModulationsNormalized(stimulusContrasts<0) = stimulusContrasts(stimulusContrasts<0)./decNormFactor;

% Compute stimulus angles in 2-D space. Let's keep angles greater than or
% equal to 0 and less than 360, just to have a clear convention. 
stimAnglesRaw = zeros(length(stimulusContrasts),1);
for n = 1:length(stimAnglesRaw)
    % Figure out whether each stimulus is increment, decrement, or blank.
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
    
    % Get angle.  Set blank as NaN 
    stimAnglesRaw(n) = round(atan2d(stimulusModulationsNormalized(n,2),stimulusModulationsNormalized(n,1)));

    if (strcmp(stim1,'Blank') && strcmp(stim2,'Blank'))
        stimAnglesRaw(n) = NaN;
    end
end
while (any(stimAnglesRaw < 0))
    stimAnglesRaw(stimAnglesRaw < 0) = stimAnglesRaw(stimAnglesRaw < 0) + 360;
end
while (any(stimAnglesRaw >= 360))
    stimAnglesRaw(stimAnglesRaw >= 360) = stimAnglesRaw(stimAnglesRaw >= 360) - 360;
end

%% Reflect around 45 degree line if desired
stimAngles = stimAnglesRaw;
if (options.refl)
    for n = 1:length(stimAnglesRaw)
        if (stimAnglesRaw(n) > 45 && stimAnglesRaw(n) < 225)
            delta = stimAnglesRaw(n) - 45;
            stimAngles(n) = 45 - delta;
        end
    end
end
while (any(stimAngles < 0))
    stimAngles(stimAngles < 0) = stimAngles(stimAngles < 0) + 360;
end
while (any(stimAngles >= 360))
    stimAngles(stimAngles >= 360) = stimAngles(stimAngles >= 360) - 360;
end
% figure; 
% plot(stimAnglesRaw,stimAngles,'ro','MarkerFaceColor','r','MarkerSize',8);

% Fix up angles for 9/14/21 small size data
% for n = 1:length(stimAngles)
%     if (stimAngles(n) >= 329 & stimAngles(n) <= 338)
%         stimAngles(n) = 334;
%     end
% end

%% Compute the false positive rate
%
% Values of 3 and 1 are from the gamepad controller used in the experiment
yesVal = 3;
noVal = 1;
numCatchTrials = length(responseVector(mean(stimulusModulationsNormalized,2)==0));
catchResponses = responseVector(mean(stimulusModulationsNormalized,2)==0);
numCatchPos = length(catchResponses(catchResponses==yesVal));
fprintf('False positive rate: %d out of %d; (%.2f percent)\n', numCatchPos, numCatchTrials, 100*(numCatchPos/numCatchTrials));

%% Compute PF for each modulation direction individually
%
% Convert absolute modulations into vector distance from the origin; can be
% thought of ans an absolute modulation "magnitude"
fprintf('Fitting modulation directions individually... ');
stimAngleList = unique(stimAngles(isfinite(stimAngles)));
modulationVectorLengths = sqrt(stimulusModulationsNormalized(:,1).^2 + stimulusModulationsNormalized(:,2).^2);

% Pre-allocate
paramsFitted_Individual = nan(length(stimAngleList),4);

% For PF evaluation
xEval = linspace(0,10.^0.5,1000);

% Plot into this figure
figure; hold on;
set(gcf, 'Color', 'w', 'Units', 'inches', 'Position', [1 1 12 4]);

% Go through all angles
for angleNum = 1:length(stimAngleList)
    % Get data for this angle
    fitIndex = find(stimAngles==stimAngleList(angleNum));
    fitResponses = responseVector(fitIndex);
    fitModulationVectorLengths = modulationVectorLengths(fitIndex);

%     % Fix up levels for 9/14/21 small size data
%     temp = unique(fitModulationVectorLengths);
%     if (length(find(fitModulationVectorLengths == temp(end))) == 20)
%         index = find(fitModulationVectorLengths == temp(end));
%         for jj = 1:length(index)/2
%             fitModulationVectorLengths(index(jj)) = 0.999*fitModulationVectorLengths(index(jj));
%         end
%     end

    if (angleNum == 1)
        numberUniqueVectorLengths = length(unique(fitModulationVectorLengths));
    else
        if (length(unique(fitModulationVectorLengths)) ~= numberUniqueVectorLengths)
            error('Assumption of same number of stimulus levels in each director violated');
        end
    end
    
    % Group data (not necessary but makes for easier plotting)
    % 
    % Inserting 0.001 as value for catch trials
    stimLevels(angleNum,:) = [0.001 unique(fitModulationVectorLengths)']; %#ok<SAGROW>
    outOfNum(angleNum,:) = zeros(1,size(stimLevels,2)); %#ok<SAGROW>
    numPos(angleNum,:) = zeros(1,size(stimLevels,2)); %#ok<SAGROW>
    
    % Build up data for each fit
    for j = 1:length(stimLevels(angleNum,:))
        if j == 1 % Catch trial
            outOfNum(angleNum,j) = numCatchTrials;
            numPos(angleNum,j) = numCatchPos;
        else
            tempRespVector = fitResponses(fitModulationVectorLengths==stimLevels(angleNum,j));
            outOfNum(angleNum, j) = length(tempRespVector);
            numPos(angleNum, j) = length(tempRespVector(tempRespVector==yesVal));
        end
    end
    
    if (options.corrGuess)
        pFA = numCatchPos/numCatchTrials;
        pHit(angleNum,:) = numPos(angleNum,:)./outOfNum(angleNum,:);
        pFACorrect = CorrectForGuessing(pFA,pFA);
        pHitCorrected(angleNum,:) = CorrectForGuessing(pHit(angleNum,:),pFA);
        index = find(pHitCorrected(angleNum,:) < 0);
        pHitCorrected(angleNum,index) = 0;
        numCatchPosFit = round(pFACorrect*numCatchTrials);
        numPosFit(angleNum,:) = round(pHitCorrected(angleNum,:) .* outOfNum(angleNum,:));
    else
        numCatchPosFit = numCatchPos;
        numPosFit = numPos;
    end
    
    % Set up the PF fitting (requires Palamedes toolbox)
    PF = @PAL_Logistic;
    falsePosProp = numCatchPosFit/numCatchTrials;
    searchGrid = [log10(mean(modulationVectorLengths)) 5 falsePosProp, 0.01];
    paramsFree = [1 1 0 0]; %[thresh slope guess lapse]; 0 = fixed; 1 = free
    paramsFitted_Individual(angleNum,:) = PAL_PFML_Fit([-3 log10(stimLevels(angleNum,:))], [numCatchPosFit numPosFit(angleNum,:)], [numCatchTrials outOfNum(angleNum,:)], searchGrid, paramsFree, PF);
    
    % Add to plot
    ax1 = subplot(1,3,1); hold on;
    h = plot(log10(xEval), PF(paramsFitted_Individual(angleNum,:), log10(xEval)), '-', 'LineWidth', 2);
    plot(log10(stimLevels(angleNum,:)), numPosFit(angleNum,:)./outOfNum(angleNum,:), 's', ...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor', h.Color, 'MarkerSize', 10, 'HandleVisibility', 'off');
end

% Finish up plot
legend(num2str(stimAngleList), 'Location', 'NorthWest')
xlim([-1 0.5])
ylim([0 1]);
axis square
xlabel('Log10 modulation intensity (au)', 'FontSize', 14);
ylabel('Prop seen', 'FontSize', 14);
title('Slopes unconstrained');
fprintf('Done.\n');

%% Fit everything together with slopes and guess/lapse rates to be equal across stimulus angles
fprintf('Now fitting all angles with the same slope and guess rate... ');
paramsFitted_Multi = PAL_PFML_FitMultiple(log10(stimLevels), numPosFit, outOfNum, searchGrid, PF, 'slopes', 'constrained', 'guessrates', ...
    'fixed', 'lapserates', 'constrained', 'lapselimits', [0 0.05]);
fprintf('Done.\n');

% Plot
for angleNum = 1:length(stimAngleList)
    ax2 = subplot(1,3,2); hold on;
    h2 = plot(log10(xEval), PF(paramsFitted_Multi(angleNum,:), log10(xEval)), '-', 'LineWidth', 2);
    subplot(1,3,2), plot(log10(stimLevels(angleNum,:)), numPosFit(angleNum,:)./outOfNum(angleNum,:), 's', 'MarkerFaceColor', h2.Color, ...
        'MarkerEdgeColor', 'none', 'MarkerSize', 10, 'HandleVisibility', 'off');
end
legend(num2str(stimAngleList), 'Location', 'NorthWest')
xlim([-1 0.5])
ylim([0 1]);
title('Slopes constrained');
axis square;
xlabel('Log10 modulation intensity (au)', 'FontSize', 14);
ylabel('Prop seen', 'FontSize', 14);

%% Evaluate fit a specific prop seen and plot on 2-D modulation space
ax3 = subplot(1,3,3); hold on;
plot([0 0], [-1.5 1.5], 'k-', 'LineWidth', 1.5);
plot( [-1.5 1.5],[0 0], 'k-', 'LineWidth', 1.5);

%% Evaluate the PF and convert back to x,y coordinates for ellipse fitting
%
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
axLim = 1.5;
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
save(fullfile(analysisDir,sprintf('%s_incDecFits_ConstrainedSlope.mat', options.subj)), 'stimAngleList', 'falsePosProp', 'paramsFitted_Individual', 'paramsFitted_Multi', 'PF');
print(gcf, fullfile(analysisDir,sprintf('%s_incDecFits_ConstrainedSlope.tiff', options.subj)), '-dtiff');

end

%% Correct for guessing
function [pCorrected] = CorrectForGuessing(pHit,pFA)

pCorrected = (pHit-pFA)/(1-pFA);

end
