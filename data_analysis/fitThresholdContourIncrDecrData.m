function fitThresholdContourIncrDecrData(options)
% fitThresholdContourIncrDecrData
%
% Script to combine data from inc/dec experiments and fit PFs to various
% axes in the 2D stimulus space
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
%      'reflIn'     - Boolean. Read fits that treat stim1 and stim2 as
%                     symmetric. Default false.
%      'reflOut'    - Boolean. Reflect data before analysis. Default false.
%      'constraintedSlopeFits' - Boolean. Fit data with PF constrained
%                     slopes. Default true
%      'scaleDecr'  - Boolean. Scale decrement thresholds to match
%                     increments. Default false.
%      'lockSlope'  - Boolean. Lock the slope of the linear fit to -1.
%                     Default true.
%
% See also: FitEllipseQ, PointsOnEllipseQ, EllipsoidMatricesGenerate.
%

% History
%   03/07/21  dhb  Cleaning up.

%% Pick up optional arguments
arguments
    options.subj string = '11043';
    options.dataDate string = '20200131';
    options.subProject string = 'IncrDecr1';
    options.condition string = 'Separation_1';
    options.norm (1,1) logical = false;
    options.corrGuess (1,1) logical = true;
    options.reflIn (1,1) logical = false;
    options.reflOut (1,1) logical = false;
    options.constrainedSlopeFits (1,1) logical = true;
    options.scaleDecr (1,1) logical = false;
    options.lockSlope (1,1) logical = true;
end

%% Housekeeping
close all
baseProject = 'AOPsychophysics';
subProject = options.subProject;

%% Set input/output name strings
if (options.scaleDecr)
    scaleDecrStr = 'scaleDecr';
else
    scaleDecrStr = 'noScaleDecr';
end
if (options.norm)
    normStr = 'norm';
else
    normStr = 'notnorm';
end
if (options.corrGuess)
    corrGuessStr = 'corrguess';
else
    corrGuessStr = 'notcorrguess';
end
if (options.reflIn)
    reflInStr = 'refl';
else
    reflInStr = 'noRefl';
end
if (options.reflOut)
    reflOutStr = 'reflOut';
else
    reflOutStr = 'noReflOut';
end

%% If fit gets stuck, try futzing with the value of errorScalar.
errorScalar = 1000;

%% Which fits?
%
% To normalize, or not to normalize? Select false to work with the fits
% from the raw modulations, true for those from the normalized ones.
if (options.norm)
    error('Need to implement analysis of fits from normalized modulations');
end

% Analyze fits with slopes constrained?
if (~options.constrainedSlopeFits)
    error('Need to implement analysis of unconstrained slope fits');
end

%% Load in the data files (update directories to wherever these data live on your machine)
%
% Select subject
%   subj = '11043'; % WST
%   subj = '11046'; % DHB
switch (options.subj)
    case {'11043' '11046'}
        fprintf('Subject ID: %s\n', options.subj);
    otherwise
        error('Specified subject number invalid')
end

%% Set up directories.
%
% Note that we're reading the output of the program
% that fits psychometric functions, so the data comes from 'analysisDir'.
analysisSubDir = sprintf('%s_%s_%s',normStr,corrGuessStr,reflInStr);
analysisBaseDir = getpref(baseProject,'analysisDir');
analysisDir = fullfile(analysisBaseDir,subProject,options.subj,options.dataDate,options.condition,analysisSubDir);
analysisOutDir = fullfile(analysisDir,sprintf('%s_%s',scaleDecrStr,reflOutStr));
if (~exist(analysisOutDir,'dir'))
    mkdir(analysisOutDir);
end

%% Read output of psychometric fitting
theData = load(fullfile(analysisDir,sprintf('%s_incDecFits_ConstrainedSlope.mat',options.subj)));
if isempty(theData)
    error('No fit data found');
end

%% Evaluate the PF and convert back to x,y coordinates for ellipse fitting
%
% First, round to the nearest .1 above the guess rate (code won't be happy if you
% try to evaluate the PF for prop seen below the lower asymptote)
PF = @PAL_Logistic;
propEvalStart = ceil(10*theData.falsePosProp)/10;
propSeen_Fit = propEvalStart:.1:.9;
modLevels_PF = nan(size(theData.paramsFitted_Multi,1), length(propSeen_Fit));
for n = 1:size(theData.paramsFitted_Multi,1)
    modLevels_PF(n,:) = 10.^PF(theData.paramsFitted_Multi(n,:), propSeen_Fit, 'inv');
end

%% Get increment and decrement thresholds
incrThresh = [];
decrThresh = [];
for pp = 1:length(theData.stimAngleList)
    if (theData.stimAngleList(pp) == 0)
        incrThresh = [incrThresh abs(cosd(theData.stimAngleList(pp)).*modLevels_PF(pp,:))];
    elseif (theData.stimAngleList(pp) == 90)
        incrThresh = [incrThresh abs(sind(theData.stimAngleList(pp)).*modLevels_PF(pp,:))];
    elseif (theData.stimAngleList(pp) == 180)
        decrThresh = [decrThresh abs(cosd(theData.stimAngleList(pp)).*modLevels_PF(pp,:))];
    elseif (theData.stimAngleList(pp) == 270)
        decrThresh = [decrThresh abs(sind(theData.stimAngleList(pp)).*modLevels_PF(pp,:))];
    end
end
incrThresh = mean(incrThresh);
decrThresh = mean(decrThresh);

%% Convert to x and y coordinates
%
% As we do this, we also reflect the symmetric data points,
% which makes the plots look better and also lets us
% make some didactically useful plots.
%
% The angles are nominal if you scale the decrements to match the
% increments, but they are only used to pull out points for explanatory
% plots so this is OK.
index = 1;
for pp = 1:length(theData.stimAngleList)
    % Get data in x,y form
    stimXDataFit(index,:) = cosd(theData.stimAngleList(pp)).*modLevels_PF(pp,:);
    stimYDataFit(index,:) = sind(theData.stimAngleList(pp)).*modLevels_PF(pp,:);
    stimAnglesFit(index) = theData.stimAngleList(pp);
    index = index+1;
    
    % This does the reflection.  Double up points on diagonal so that each
    % datum counts same number of times in the fit.
    if (options.reflOut)
        stimXDataFit(index,:) = stimYDataFit(index-1,:);
        stimYDataFit(index,:) = stimXDataFit(index-1,:);
        if (stimXDataFit(index,end) == 0 & stimYDataFit(index,end) == 0)
            error('Angle undefined the way we are computing it. Rethink reflOut logic.')
        end
        stimAnglesFit(index) = atan2d(stimYDataFit(index,end),stimXDataFit(index,end));
        index = index+1;
    end
end
while (any(stimAnglesFit < 0))
    stimAnglesFit(stimAnglesFit < 0) = stimAnglesFit(stimAnglesFit < 0) + 360;
end
while (any(stimAnglesFit >= 360))
    stimAnglesFit(stimAnglesFit >= 360) = stimAnglesFit(stimAnglesFit >= 360) - 360;
end

%% Make increments and decrements same scale, if desired
if (options.scaleDecr)
    incrDecrScaleFactor = incrThresh./decrThresh;
    for pp = 1:size(stimXDataFit,1)
        for jj = 1:size(stimXDataFit,2)
            if (stimXDataFit(pp,jj) < 0)
                stimXDataFit(pp,jj) = incrDecrScaleFactor*stimXDataFit(pp,jj);
            end
            if (stimYDataFit(pp,jj) < 0)
                stimYDataFit(pp,jj) = incrDecrScaleFactor*stimYDataFit(pp,jj);
            end
        end
    end
end

%% Get data for desired proportions correct right format.
%
% If there is more than one propsSeen specified, data go into a cell array.
% Otherwise just a regular array. The code should work even with a cell
% array of length 1, but doing it like this lets us test FitEllipseQ a
% little more thorougly.
propsSeen = 0.70;
whichCol = find(propSeen_Fit == propsSeen);
xData = stimXDataFit(:,whichCol)';
yData = stimYDataFit(:,whichCol)';
theDataToFit = [xData ; yData];

%% Plot the points
theColors = ['r' 'k' 'b' 'b' 'y' 'c'];
if (length(propsSeen) > length(theColors))
    error('Need to specify more plot colors, fewer propsSeen, or fix code to handle');
end
theEllipseFig = figure; clf; hold on;
plot(theDataToFit(1,:),theDataToFit(2,:),[theColors(1) 'o'],'MarkerFaceColor',theColors(1),'MarkerSize',12);
title(sprintf('Subject %s, Criterion %d%%',options.subj,round(100*propsSeen)));
theLim = 2;
plot([-theLim theLim],[0 0],'k:','LineWidth',1);
plot([0 0],[-theLim theLim],'k:','LineWidth',1);
xlim([-theLim theLim]);
ylim([-theLim theLim]);
axis('square');
xlabel('Contrast 1')
ylabel('Contrast 2');
print(theEllipseFig, fullfile(analysisOutDir,sprintf('%s_AllData.tiff', options.subj)), '-dtiff');

%% Fit ellipse
%
% The work of generating the points to plot is done by PointsOnEllipseQ.
[ellParams,A,Ainv,Q] = FitEllipseQ(theDataToFit,'lockAngleAt0',false,'errorScalar',errorScalar,'initialParams',[incrThresh(end) decrThresh(end) 45]);
nCirclePoints = 100;
circlePoints = UnitCircleGenerate(nCirclePoints);
ellPoints = PointsOnEllipseQ(Q,circlePoints);

% Print ellipse parameters
%
% Give full length of axes, not radius equivalent (which is why it is 2 in
% the numerator of the numbers below, not 1).  Printout is a little
% different depending on whether angle is locked or not.
fprintf('Major axis length: %0.2f, minor axis length: %0.2f, major axis angle %0.1f deg\n',2/ellParams(1),2/ellParams(2),ellParams(3));

% Add ellipse to plot
theEllipseFig1 = theEllipseFig.copy;
figure(theEllipseFig);
plot(ellPoints(1,:),ellPoints(2,:),'r','LineWidth',2);
title(sprintf('Subject %s, Criterion %d%%',options.subj,round(100*propsSeen)));
theLim = 2;
plot([-theLim theLim],[0 0],'k:','LineWidth',1);
plot([0 0],[-theLim theLim],'k:','LineWidth',1);
xlim([-theLim theLim]);
ylim([-theLim theLim]);
axis('square');
xlabel('Contrast 1')
ylabel('Contrast 2');
print(theEllipseFig, fullfile(analysisOutDir,sprintf('%s_EllipseAllData.tiff', options.subj)), '-dtiff');

%% Fit ellipse with angle locked to just pure increment and decrement data
neededAngles = [0 90 180 270];
index = [];
for aa = 1:length(neededAngles)
    index = [index find(stimAnglesFit == neededAngles(aa))];
end
[ellParams0,A0,Ainv0,Q0] = FitEllipseQ(theDataToFit(:,index),'lockAngleAt0',true,'errorScalar',errorScalar,'initialParams',[incrThresh(end) decrThresh(end) 0]);
ellPoints0 = PointsOnEllipseQ(Q0,circlePoints);

% Plot
theEllipseFig0 = figure; clf; hold on;
plot(theDataToFit(1,index),theDataToFit(2,index),[theColors(1) 'o'],'MarkerFaceColor',theColors(1),'MarkerSize',12);
title(sprintf('Subject %s, Criterion %d%%',options.subj,round(100*propsSeen)));
theLim = 2;
plot([-theLim theLim],[0 0],'k:','LineWidth',1);
plot([0 0],[-theLim theLim],'k:','LineWidth',1);
xlim([-theLim theLim]);
ylim([-theLim theLim]);
axis('square');
xlabel('Contrast 1')
ylabel('Contrast 2');
print(theEllipseFig0, fullfile(analysisOutDir,sprintf('%s_IncrAndDecrOnlyData.tiff', options.subj)), '-dtiff');
plot(ellPoints0(1,:),ellPoints0(2,:),'r','LineWidth',2);
print(theEllipseFig0, fullfile(analysisOutDir,sprintf('%s_Ellipse0FitData.tiff', options.subj)), '-dtiff');

% Add to plot of all data
figure(theEllipseFig1);
plot(ellPoints0(1,:),ellPoints0(2,:),'r','LineWidth',2);
title(sprintf('Subject %s, Criterion %d%%',options.subj,round(100*propsSeen)));
theLim = 2;
plot([-theLim theLim],[0 0],'k:','LineWidth',1);
plot([0 0],[-theLim theLim],'k:','LineWidth',1);
xlim([-theLim theLim]);
ylim([-theLim theLim]);
axis('square');
xlabel('Contrast 1')
ylabel('Contrast 2');
print(theEllipseFig1, fullfile(analysisOutDir,sprintf('%s_Ellipse0AllData.tiff', options.subj)), '-dtiff');

%% Build up some explanatory plots
%
% Just incremental stimulus one
theIncrFig = figure; clf; hold on
index = find(stimAnglesFit == 0);
if (iscell(theDataToFit))
    for pp = 1:length(propsSeen)
        plot(theDataToFit{pp}(1,index),theDataToFit{pp}(2,index),[theColors(pp) 'o'],'MarkerFaceColor',theColors(pp),'MarkerSize',12);
    end
else
    plot(theDataToFit(1,index),theDataToFit(2,index),[theColors(1) 'o'],'MarkerFaceColor',theColors(1),'MarkerSize',12);
end
plot([-theLim theLim],[0 0],'k:','LineWidth',1);
plot([0 0],[-theLim theLim],'k:','LineWidth',1);
xlim([-theLim theLim]);
ylim([-theLim theLim]);
axis('square');
xlabel('Contrast 1')
ylabel('Contrast 2');
title(sprintf('Subject %s, Criterion %d%%',options.subj,round(100*propsSeen)));
print(theIncrFig, fullfile(analysisOutDir,sprintf('%s_SingleIncr.tiff', options.subj)), '-dtiff');

%% Both incremental stimuli
index = find(stimAnglesFit == 90);
if (iscell(theDataToFit))
    for pp = 1:length(propsSeen)
        plot(theDataToFit{pp}(1,index),theDataToFit{pp}(2,index),[theColors(pp) 'o'],'MarkerFaceColor',theColors(pp),'MarkerSize',12);
    end
else
    plot(theDataToFit(1,index),theDataToFit(2,index),[theColors(1) 'o'],'MarkerFaceColor',theColors(1),'MarkerSize',12);
end
print(theIncrFig, fullfile(analysisOutDir,sprintf('%s_DoubleIncr.tiff', options.subj)), '-dtiff');
theIncrFig1 = theIncrFig.copy;
figure(theIncrFig);

%% Fit linear model to single increments and plot
neededAngles = [0 90];
index = [];
for aa = 1:length(neededAngles)
    index = [index find(stimAnglesFit == neededAngles(aa))];
end
x = theDataToFit(1,index)';
y = theDataToFit(2,index)';
lineParams = [x ones(size(x))]\y;
linePlotX = linspace(-theLim,theLim,100)';
linePlotY = [linePlotX ones(size(linePlotX))]*lineParams;
plot(linePlotX,linePlotY,'r','LineWidth',3);
print(theIncrFig, fullfile(analysisOutDir,sprintf('%s_DoubleIncrWithLine.tiff', options.subj)), '-dtiff');

%% Add incr-incr datum
index = find(stimAnglesFit == 45);
plot(theDataToFit(1,index),theDataToFit(2,index),[theColors(1) 'o'],'MarkerFaceColor',theColors(1),'MarkerSize',12);
print(theIncrFig, fullfile(analysisOutDir,sprintf('%s_IncrWithDoubleIncrLine.tiff', options.subj)), '-dtiff');

%% Add incr-incr
figure(theIncrFig1);
index = find(stimAnglesFit == 45);
plot(theDataToFit(1,index),theDataToFit(2,index),[theColors(1) 'o'],'MarkerFaceColor',theColors(1),'MarkerSize',12);
print(theIncrFig1, fullfile(analysisOutDir,sprintf('%s_Incr.tiff', options.subj)), '-dtiff');

%% Fit linear model to all the incremental data and plot
neededAngles = [0 45 90];
index = [];
for aa = 1:length(neededAngles)
    index = [index find(stimAnglesFit == neededAngles(aa))];
end
x = theDataToFit(1,index)';
y = theDataToFit(2,index)';
if (options.lockSlope)
    slope = -1;
    lineParams = [slope ones(size(x))\(y-(slope*x))]';
else
    lineParams = [x ones(size(x))]\y;
end
linePlotX = linspace(-theLim,theLim,100)';
linePlotY = [linePlotX ones(size(linePlotX))]*lineParams;
plot(linePlotX,linePlotY,'r','LineWidth',3);
print(theIncrFig1, fullfile(analysisOutDir,sprintf('%s_IncrWithBestIncrLine.tiff', options.subj)), '-dtiff');

%% Add decr-decr points
neededAngles = [180 225 270];
index = [];
for aa = 1:length(neededAngles)
    index = [index find(stimAnglesFit == neededAngles(aa))];
end
plot(theDataToFit(1,index),theDataToFit(2,index),[theColors(1) 'o'],'MarkerFaceColor',theColors(1),'MarkerSize',12);
print(theIncrFig1, fullfile(analysisOutDir,sprintf('%s_IncrAndDecrWithBestIncrLine.tiff', options.subj)), '-dtiff');

%% Add decr-decr line
neededAngles = [180 225 270];
index = [];
for aa = 1:length(neededAngles)
    index = [index find(stimAnglesFit == neededAngles(aa))];
end
x = theDataToFit(1,index)';
y = theDataToFit(2,index)';
if (options.lockSlope)
    slope = -1;
    lineParams = [slope ones(size(x))\(y-(slope*x))]';
else
    lineParams = [x ones(size(x))]\y;
end
linePlotX = linspace(-theLim,theLim,100)';
linePlotY = [linePlotX ones(size(linePlotX))]*lineParams;
plot(linePlotX,linePlotY,'r','LineWidth',3);
print(theIncrFig1, fullfile(analysisOutDir,sprintf('%s_IncrAndDecrWithBestIncrAndDecrLines.tiff', options.subj)), '-dtiff');

%% Add the critical incr-decr points
index = find((stimAnglesFit > 90 & stimAnglesFit < 180) | (stimAnglesFit > 270 & stimAnglesFit < 360));
plot(theDataToFit(1,index),theDataToFit(2,index),[theColors(1) 'o'],'MarkerFaceColor',theColors(1),'MarkerSize',12);
print(theIncrFig1, fullfile(analysisOutDir,sprintf('%s_AllDataWithBestIncrAndDecrLines.tiff', options.subj)), '-dtiff');
