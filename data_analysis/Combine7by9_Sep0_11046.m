% Combine7by9_Sep0_11046
%
% High level script to combine 7x9 data across sessions
%
% This version for the 7x9 pixel spots at 0 separation, for
% subject 11046.

%% Clear
clear; close all;

%% Specify sessions
whichSize = '7x9';
switch (whichSize)
    case '7x9'
        theFiles{1} = '/Users/dhb/Dropbox (Aguirre-Brainard Lab)/AOPY_analysis/AOPsychophysics/IncrDecr1/11046/20200131/Separation_1/notnorm_corrguess_norefl/scaleDecr_noReflOut/11046_ContourAnalysis.mat';
        theFiles{2} = '/Users/dhb/Dropbox (Aguirre-Brainard Lab)/AOPY_analysis/AOPsychophysics/IncrDecr2/11046/20210914/Size_1/notnorm_corrguess_noRefl/scaleDecr_noReflOut/11046_ContourAnalysis.mat';
        theFiles{3} = '/Users/dhb/Dropbox (Aguirre-Brainard Lab)/AOPY_analysis/AOPsychophysics/IncrDecr4/11046/20211026/Size_1/notnorm_corrguess_norefl/scaleDecr_noReflOut/11046_ContourAnalysis.mat';
        theFiles{4} = '/Users/dhb/Dropbox (Aguirre-Brainard Lab)/AOPY_analysis/AOPsychophysics/IncrDecr5/11046/20211123/Size1_Sep0/notnorm_corrguess_norefl/noScaleDecr_noReflOut/11046_ContourAnalysis.mat';
        titleStr = '11046, 7x9';
        theLim = 2;
end

%% Load in data from each session
theDataToFit = [];
stimAnglesFit = [];
for ii = 1:length(theFiles)
    theData{ii} = load(theFiles{ii},'theDataToFit','stimAnglesFit');
    theDataToFit = [theDataToFit theData{ii}.theDataToFit];
    stimAnglesFit = [stimAnglesFit theData{ii}.stimAnglesFit];
end

%% Plot the data
theColors = ['r' 'g' 'b' 'k' 'c' 'y'];
theColor = 1;
if (length(theData) > length(theColors))
    error('Need to specify more plot colors, fewer propsSeen, or fix code to handle');
end
theDataFig = figure; clf; hold on;
for ii = 1:length(theData)
    plot(theData{ii}.theDataToFit(1,:),theData{ii}.theDataToFit(2,:),[theColors(ii) 'o'],'MarkerFaceColor',theColors(ii),'MarkerSize',12);

    % Cycle through colors
    theColor = theColor + 1;
    if (theColor > length(theColors))
        theColor = 1;
    end
end
title([titleStr ' Not Scaled']);
plot([-theLim theLim],[0 0],'k:','LineWidth',1);
plot([0 0],[-theLim theLim],'k:','LineWidth',1);
xlim([-theLim theLim]);
ylim([-theLim theLim]);
axis('square');
xlabel('Contrast 1')
ylabel('Contrast 2');
% print(theDataFig, fullfile(analysisOutDir,sprintf('%s_AllData.tiff', options.subj)), '-dtiff');

%% Read in ideal observer threshold contour and add to plot, after scaling
baseProject = 'AOCompObserver';
compAnalysisBaseDir = getpref(baseProject,'analysisDir');
computationalName = '7_9_0';
defocusDiopters = 0.10;
pupilDiam = 7;
compAnalysisInDir = fullfile(compAnalysisBaseDir,sprintf('%s_%s_%d',computationalName,num2str(round(1000*defocusDiopters)),pupilDiam));
if (~exist(compAnalysisInDir ,'dir'))
    error('Computational observer not yet run for specified diopters of defocus');
end
compObserver = load(fullfile(compAnalysisInDir,sprintf('CompObserver_%s',computationalName)));

nCirclePoints = 100;
circlePoints = UnitCircleGenerate(nCirclePoints);

circlePointsFit(1,:) = cosd(stimAnglesFit);
circlePointsFit(2,:) = sind(stimAnglesFit);
compObserverEllData = PointsOnEllipseQ(compObserver.compFitQ,circlePointsFit);
for aa = 1:length(stimAnglesFit)
    dataRadii(aa) = norm(theDataToFit(:,aa));
    compRadii(aa) = norm(compObserverEllData(:,aa));
end
compFitFactor = compRadii'\dataRadii';
compObserverEll = PointsOnEllipseQ(compObserver.compFitQ,circlePoints)*compFitFactor;
%theEllipseFig3 = figure; hold on
plot(compObserverEll(1,:),compObserverEll(2,:),'r','LineWidth',2);
%title(sprintf('Subject %s, Criterion %d%%, Fit w/ Ideal, %0.2g D',options.subj,round(100*propsSeen),options.defocusDiopters));
fprintf('Ideal observer scale factor: %0.3g\n',compFitFactor);
plot([-theLim theLim],[0 0],'k:','LineWidth',1);
plot([0 0],[-theLim theLim],'k:','LineWidth',1);
xlim([-theLim theLim]);
ylim([-theLim theLim]);
axis('square');
xlabel('Contrast 1')
ylabel('Contrast 2');
%print(theEllipseFig3, fullfile(analysisOutDir,sprintf('%s_%s_CompEllipse.tiff',options.subj,num2str(round(1000*options.defocusDiopters)))), '-dtiff');
%print(theEllipseFig3, fullfile(analysisOutDir,sprintf('%s_%s_CompEllipseWithData.tiff',options.subj,num2str(round(1000*options.defocusDiopters)))), '-dtiff');

% %% Plot data scaled to first session 
% %
% % This code only works because either a session has data on all of 0, 90, 180 and
% % 270, or on none of these.  A more general bit of code would need to find
% % the common angles across sessions in some way and use those to scale.
% %
% % Probably scaling isn't really needed and might not be a good idea anyway.
% theColors = ['r' 'g' 'b' 'k' 'c' 'y'];
% theColor = 1;
% theDataScaledFig = figure; clf; hold on;
% for ii = 1:length(theData)
%     % Don't scale first session, this is serving as a reference
%     if (ii == 1)
%         xScaleIndex = find(theData{ii}.theDataToFit(2,:) == 0);
%         yScaleIndex = find(theData{ii}.theDataToFit(1,:) == 0);
%         scaleData(:,ii) = [theData{ii}.theDataToFit(1,xScaleIndex) theData{ii}.theDataToFit(2,yScaleIndex)]';
%         scaleFactor(ii) = 1;
%     else
%         xScaleIndex = find(theData{ii}.theDataToFit(2,:) == 0);
%         yScaleIndex = find(theData{ii}.theDataToFit(1,:) == 0);
%         
%         % Can only scale if there was on-axis data for the session.  This
%         % is not always true.
%         if isempty(xScaleIndex) & isempty(yScaleIndex)
%             scaleFactor(ii) = 1;
%         else
%             scaleData(:,ii) = [theData{ii}.theDataToFit(1,xScaleIndex) theData{ii}.theDataToFit(2,yScaleIndex)]';
%             scaleFactor(ii) = scaleData(:,ii)\scaleData(:,1);
%         end
%     end
% 
%     % Add to plot
%     plot(scaleFactor(ii)*theData{ii}.theDataToFit(1,:),scaleFactor(ii)*theData{ii}.theDataToFit(2,:),[theColors(ii) 'o'],'MarkerFaceColor',theColors(ii),'MarkerSize',12);
% 
%     % Cycle through colors
%     theColor = theColor + 1;
%     if (theColor > length(theColors))
%         theColor = 1;
%     end
% end
% title([titleStr ' Scaled']);
% plot([-theLim theLim],[0 0],'k:','LineWidth',1);
% plot([0 0],[-theLim theLim],'k:','LineWidth',1);
% xlim([-theLim theLim]);
% ylim([-theLim theLim]);
% axis('square');
% xlabel('Contrast 1')
% ylabel('Contrast 2');
% % print(theDataFig, fullfile(analysisOutDir,sprintf('%s_AllData.tiff', options.subj)), '-dtiff');
