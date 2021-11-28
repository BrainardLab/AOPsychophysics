% Combine7by9
%
% High level script to combine 7x9 data across sessions

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
%  
end

%% Load in data from each session
for ii = 1:length(theFiles)
    theData{ii} = load(theFiles{ii},'theDataToFit');
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

%% Plot data scaled to first session 
%
% The fourth session just has one direction, so we can't scale it
theColors = ['r' 'g' 'b' 'k' 'c' 'y'];
theColor = 1;
theDataScaledFig = figure; clf; hold on;
for ii = 1:length(theData)
    
    if (ii == 1)
        xScaleIndex = find(theData{ii}.theDataToFit(2,:) == 0);
        yScaleIndex = find(theData{ii}.theDataToFit(1,:) == 0);
        scaleData(:,ii) = [theData{ii}.theDataToFit(1,xScaleIndex) theData{ii}.theDataToFit(2,yScaleIndex)]';
        scaleFactor(ii) = 1;
        plot(theData{ii}.theDataToFit(1,:),theData{ii}.theDataToFit(2,:),[theColors(ii) 'o'],'MarkerFaceColor',theColors(ii),'MarkerSize',12);
    elseif (ii == 4)
        scaleFactor(ii) = 1;
        plot(scaleFactor(ii)*theData{ii}.theDataToFit(1,:),scaleFactor(ii)*theData{ii}.theDataToFit(2,:),[theColors(ii) 'o'],'MarkerFaceColor',theColors(ii),'MarkerSize',12); 
    else
        xScaleIndex = find(theData{ii}.theDataToFit(2,:) == 0);
        yScaleIndex = find(theData{ii}.theDataToFit(1,:) == 0);
        scaleData(:,ii) = [theData{ii}.theDataToFit(1,xScaleIndex) theData{ii}.theDataToFit(2,yScaleIndex)]';
        scaleFactor(ii) = scaleData(:,ii)\scaleData(:,1);
        plot(scaleFactor(ii)*theData{ii}.theDataToFit(1,:),scaleFactor(ii)*theData{ii}.theDataToFit(2,:),[theColors(ii) 'o'],'MarkerFaceColor',theColors(ii),'MarkerSize',12); 
    end

    % Cycle through colors
    theColor = theColor + 1;
    if (theColor > length(theColors))
        theColor = 1;
    end
end
title([titleStr ' Scaled']);
plot([-theLim theLim],[0 0],'k:','LineWidth',1);
plot([0 0],[-theLim theLim],'k:','LineWidth',1);
xlim([-theLim theLim]);
ylim([-theLim theLim]);
axis('square');
xlabel('Contrast 1')
ylabel('Contrast 2');
% print(theDataFig, fullfile(analysisOutDir,sprintf('%s_AllData.tiff', options.subj)), '-dtiff');
