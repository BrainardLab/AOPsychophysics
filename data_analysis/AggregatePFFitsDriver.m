% AggregatePFFitsDriver
%
% Description:    
%    The initial PF fits constrain PF slope within condition, but sometimes
%    there isn't enough data to do that, and it makes sense to leverage more
%    data and get the overall fit with constrained PF slope. This script does
%    the aggregation.
%
%    Call into this after setting up analysis parameters

% History:
%    12/23/21  dhb  Pull out driver

%% Load in data from each session
warnState = warning('off','MATLAB:load:cannotInstantiateLoadedVariable');
sessionData_All = {};
stimAnglesFit = [];
dirIndices_All = [];
sessionNumbers_All = [];
checkIndex = 1;
firstSession = true;
for ii = 1:length(theFiles)
    theData{ii} = load(theFiles{ii},'sessionData');
    sessionData_All = {sessionData_All{:} theData{ii}.sessionData{:} };

    % Session number check logic
    if (firstSession)
        currentSession = sessionNames{ii};
        currentSessionNumber = sessionNumbers(ii);
        firstSession = false;
    end
    if (strcmp(sessionNames{ii},currentSession))
        if (sessionNumbers(ii) ~= currentSessionNumber)
            error('Session name/number mismatch');
        end
    else
        currentSession = sessionNames{ii};
        currentSessionNumber = sessionNumbers(ii);
    end

    % Sanity checks as well as create index from each session summary back to the
    % directory it was read from.
    for ss = 1:length(theData{ii}.sessionData)
        sessionNumbers_All = [sessionNumbers_All currentSessionNumber];
        
        % Check height
        if (theData{ii}.sessionData{ss}.stimHeight ~= stimHeightCheck(ii))
            error('Stim height check failure');
        end

        % Check width
        if (theData{ii}.sessionData{ss}.stimWidth ~= stimWidthCheck(ii))
            error('Stim width check failure');
        end

        % This check assumes only one separation per directory, which is
        % currently true of our data collection practice.
        if (theData{ii}.sessionData{ss}.stimSeparationPixels ~= stimSeparationPixelsCheck(ii))
            error('Stim separation check failure');
        end

        % Preserve indices into individual aggregated files
        dirIndices_All = [dirIndices_All ii];
    end
end
warning(warnState);

%% Get data to fit into one big matrix
for ss = 1:length(sessionData_All)
    angles_All(ss) = sessionData_All{ss}.angle;
    stimLevels_All(ss,:) = sessionData_All{ss}.stimLevels;
    outOfNum_All(ss,:) = sessionData_All{ss}.outOfNum;
    numPosFit_All(ss,:) = sessionData_All{ss}.numPosFit;
    stimHeightPixels_All(ss) = sessionData_All{ss}.stimHeight;
    stimWidthPixels_All(ss) = sessionData_All{ss}.stimWidth;
    stimSeparationPixels_All(ss) = sessionData_All{ss}.stimSeparationPixels;
end

%% Fit the whole pile of data
%
% The constraint on the guess rate at 0 makes sense given
% correction for guessing, but not more generally.  Lapse
% can be higher than usual after correction for guessing.
PF = @PAL_Logistic;
initialParams = [log10(mean(stimLevels_All(:))) 15 0 0.01];
initialAlphas = [0.5*log10(mean(stimLevels_All(:))) log10(mean(stimLevels_All(:))) 1.5*log10(mean(stimLevels_All(:)))];
initialBetas = [7.5 15 22.5];
idx = 1;
for aa = 1:length(initialAlphas)
    for bb = 1:length(initialBetas)
        initialParams(1) = initialAlphas(aa);
        initialParams(2) = initialBetas(bb);
        [paramsFittedAggregate_Temp{idx} LL(idx)] = PAL_PFML_FitMultiple(log10(stimLevels_All), numPosFit_All, outOfNum_All, initialParams, PF, 'slopes', 'constrained', 'guessrates', ...
            'fixed', 'lapserates', 'constrained', 'lapselimits', [0 0.05]);
        idx = idx+1;
    end
end
[~,whichIdx] = max(LL);
paramsFittedAggregate_All = paramsFittedAggregate_Temp{whichIdx};

%% For PF evaluation
xEval = linspace(0,10.^0.5,1000);

%% Now deal the results back out to the individual data directories
for dd = 1:length(theFiles)
    theIndices = find(dirIndices_All == dd);
    paramsFittedAggregate = paramsFittedAggregate_All(theIndices,:);
    stimAngleList = angles_All(theIndices);
    stimLevels = stimLevels_All(theIndices,:);
    outOfNum = outOfNum_All(theIndices,:);
    numPosFit = numPosFit_All(theIndices,:);
    stimHeightPixels = stimHeightPixels_All(theIndices);
    stimWidthPixels = stimWidthPixels_All(theIndices);
    stimSeparationPixels = stimSeparationPixels_All(theIndices);
    sessionNumbers = sessionNumbers_All(theIndices);
    sessionData = { sessionData_All{theIndices} };

    % Plot for this dir.  Also aggregate up info on thresholds.  Make sure
    % to clear the variables we care about before each time we run this
    % loop, because there are different numbers of angles in different
    % data files.
    pfFig = figure; hold on;
    nXSubplots = ceil(sqrt(length(stimAngleList)));
    nYSubplots = nXSubplots;
    set(gcf, 'Color', 'w', 'Units', 'inches', 'Position', [1 1 11 11]);
    whichSubplot = 1;
    clear thresholdCriteria thresholdContrasts thresholdXData thresholdYData
    for aa = 1:length(stimAngleList)
        % One PF in each subplot
        figure(pfFig);
        subplot(nXSubplots,nYSubplots,whichSubplot); hold on;
        h = plot(log10(xEval), PF(paramsFittedAggregate(aa,:), log10(xEval)), '-', 'LineWidth', 2);
        plot(log10(stimLevels(aa,:)), numPosFit(aa,:)./outOfNum(aa,:), 's', ...
            'MarkerEdgeColor', 'none', 'MarkerFaceColor', h.Color, 'MarkerSize', 10, 'HandleVisibility', 'off');
        whichSubplot = whichSubplot + 1;
        legend(num2str(stimAngleList(aa)), 'Location', 'NorthWest')
        xlim([-1 0.5])
        ylim([0 1]);
        axis square
        xlabel('Log10 modulation intensity (au)', 'FontSize', 14);
        ylabel('Prop seen', 'FontSize', 14);
        title(sprintf('Slope: %0.2f',paramsFittedAggregate(aa,2)));

        % Compute threshold data for this angle in this file
        thresholdCriteria(aa) = thresholdCriterion;
        thresholdContrasts(aa) = 10.^PF(paramsFittedAggregate(aa,:), thresholdCriteria(aa), 'inv');
        thresholdXData(aa) = thresholdContrasts(aa)*cosd(stimAngleList(aa));
        thresholdYData(aa) = thresholdContrasts(aa)*sind(stimAngleList(aa));
    end

    % Save 
    save(fullfile(sessionData_All{theIndices(1)}.dirName,sprintf('%s_%s_Aggregated.mat',theSubject,theIndicator)), ...
        'sessionData','paramsFittedAggregate','PF','stimAngleList','stimLevels','outOfNum','numPosFit','stimHeightPixels', ...
        'stimWidthPixels','stimSeparationPixels','sessionNumbers', ...
        'thresholdCriteria','thresholdContrasts','thresholdXData','thresholdYData');
    print(pfFig, fullfile(sessionData_All{theIndices(1)}.dirName,sprintf('%s_%s_Aggregated.tiff',theSubject,theIndicator)), '-dtiff');
end

%% Close
close all;

