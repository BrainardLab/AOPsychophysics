% AggregatePFFits
%
% The initial PF fits constrain PF slope within condition, but sometimes
% there isn't enough data to do that, and it makes sense to leverage more
% data and get the overall fit with constrained PF slope. This script does
% the aggregation.
%
% RunAllInitialPFFits must have been run first, to put the data in the form
% where this program can use it.

% History:
%    12/04/21  dhb  Wrote it.

%% Clear
clear; close all;

%% Parameters
%
% Where's are the initial PF fits?
psychoProject = 'AOPsychophysics';
psychoBaseDir = getpref(psychoProject,'analysisDir');

%% 11046, 7x9 data
theFiles{1} = fullfile(psychoBaseDir,'IncrDecr1','11046','20200131','Separation_1','notnorm_corrguess_norefl','11046_incDecFits_ConstrainedSlope.mat');
theFiles{2} = fullfile(psychoBaseDir,'IncrDecr2','11046','20210914','Size_1','notnorm_corrguess_norefl','11046_incDecFits_ConstrainedSlope.mat');
theFiles{3} = fullfile(psychoBaseDir,'IncrDecr4','11046','20211026','Size_1','notnorm_corrguess_norefl','11046_incDecFits_ConstrainedSlope.mat');
theFiles{4} = fullfile(psychoBaseDir,'IncrDecr5','11046','20211123','Size1_Sep0','notnorm_corrguess_norefl','11046_incDecFits_ConstrainedSlope.mat');
theFiles{5} = fullfile(psychoBaseDir,'IncrDecr5','11046','20211123','Size1_Sep2','notnorm_corrguess_norefl','11046_incDecFits_ConstrainedSlope.mat');
theFiles{6} = fullfile(psychoBaseDir,'IncrDecr5','11046','20211123','Size1_Sep4','notnorm_corrguess_norefl','11046_incDecFits_ConstrainedSlope.mat');
theFiles{7} = fullfile(psychoBaseDir,'IncrDecr5','11046','20211123','Size1_Sep6','notnorm_corrguess_norefl','11046_incDecFits_ConstrainedSlope.mat');
theFiles{8} = fullfile(psychoBaseDir,'IncrDecr5','11046','20211123','Size1_Sep8','notnorm_corrguess_norefl','11046_incDecFits_ConstrainedSlope.mat');

%% Load in data from each session
warnState = warning('off','MATLAB:load:cannotInstantiateLoadedVariable');
sessionData = {};
stimAnglesFit = [];
dirIndices = [];
for ii = 1:length(theFiles)
    theData{ii} = load(theFiles{ii},'sessionData');
    sessionData = { sessionData{:} theData{ii}.sessionData{:} };
    for ss = 1:length(theData{ii}.sessionData)
        dirIndices = [dirIndices ii];
    end
end
warning(warnState);

%% Get data to fit into one big matrix
for ss = 1:length(sessionData)
    angles_All(ss) = sessionData{ss}.angle;
    stimLevels_All(ss,:) = sessionData{ss}.stimLevels;
    outOfNum_All(ss,:) = sessionData{ss}.outOfNum;
    numPosFit_All(ss,:) = sessionData{ss}.numPosFit;
end

%% Fit the whole pile of data
%
% The constraint on the guess rate at 0 makes sense given
% correction for guessing, but not more generally.  Lapse
% can be higher than usual after correction for guessing.
PF = @PAL_Logistic;
searchGrid = [log10(mean(stimLevels_All(:))) 5 0, 0.01];
paramsFittedAggregate_All = PAL_PFML_FitMultiple(log10(stimLevels_All), numPosFit_All, outOfNum_All, searchGrid, PF, 'slopes', 'constrained', 'guessrates', ...
    'fixed', 'lapserates', 'constrained', 'lapselimits', [0 0.10]);

%% For PF evaluation
xEval = linspace(0,10.^0.5,1000);

%% Now deal the results back out to the individual data directories
for dd = 1:length(theFiles)
    theIndices = find(dirIndices == dd);
    paramsFittedAggregate = paramsFittedAggregate_All(theIndices,:);
    stimAngleList = angles_All(theIndices);
    stimLevels = stimLevels_All(theIndices,:);
    outOfNum = outOfNum_All(theIndices,:);
    numPosFit = numPosFit_All(theIndices,:);

    % Plot for this dir
    pfFig = figure; hold on;
    nXSubplots = ceil(sqrt(length(stimAngleList)));
    nYSubplots = nXSubplots;
    set(gcf, 'Color', 'w', 'Units', 'inches', 'Position', [1 1 11 11]);
    whichSubplot = 1;
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
    end

    % Save 
    save(fullfile(sessionData{theIndices(1)}.dirName,'AggregatedPFFit.mat'),'paramsFittedAggregate','PF','stimAngleList','stimLevels','outOfNum','numPosFit');
    print(pfFig, fullfile(sessionData{theIndices(1)}.dirName,'AggregatedPFPlot.tiff'), '-dtiff');
end

%% Close
close all;

