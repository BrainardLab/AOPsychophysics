function fitThresholdContourIncrDecrData(options)
% fitThresholdContourIncrDecrData
%
% Function to combine data from inc/dec experiments and fit PFs to various
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
%      'defocusDiopters' - Numeric.  Compuational observer assumed defocus
%                     in diopters. Default 0.05.
%      'computationalName' - String. Name of computational observer
%                     condition to load. Default '7_9_0'
%      'pfSlope'    - Locked slope used in PF fit.  Default empty, measn
%                     unlocked fit file used. This causes the routine to
%                     read from a file with a slope indicator in its name,
%                     and (somewhat confusingly) use the non constrained
%                     fit variable, which for this case has been done with the
%                     passed constraint.
%      'theIndictor' - String. Goes into input and output datafinames.
%                     Default 'incDecFits'
%      'aggregateFitName' - String. Use fits from aggregated data PF analysis.
%                     This is the filename to use. Default empty, in which
%                     case the default within session PF fits are used.
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
    options.defocusDiopters (1,1) double = 0.05;
    options.pupilDiam (1,1) double = 7;
    options.theLim (1,1) double = 2;
    options.computationalName = '7_9_0'
    options.pfSlope = [];
    options.theIndicator = 'incDecFits';
    options.aggregateFitName = '';
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
    reflInStr = 'norefl';
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
% that fits psychometric functions, so the data comes from 'analysisInDir'.
analysisBaseDir = getpref(baseProject,'analysisDir');
analysisSubDir = sprintf('%s_%s_%s',normStr,corrGuessStr,reflInStr);
analysisDir = fullfile(analysisBaseDir,subProject,options.subj,options.dataDate,options.condition,analysisSubDir);
if (isempty(options.aggregateFitName))
    analysisOutDir = fullfile(analysisDir,sprintf('%s_%s',scaleDecrStr,reflOutStr));
else
    analysisOutDir = fullfile(analysisDir,sprintf('%s_%s_%s',scaleDecrStr,reflOutStr,options.aggregateFitName));
end

if (~exist(analysisOutDir,'dir'))
    mkdir(analysisOutDir);
end

%% Read output of psychometric fitting
if (isempty(options.aggregateFitName))
    if (isempty(options.pfSlope))
        theData = load(fullfile(analysisDir,sprintf('%s_%s_ConstrainedSlope.mat',options.subj,options.theIndicator)));
        thePFParams = theData.paramsFitted_Multi;
    else
        % Note that the fixed slope fit parameters are in
        % theData.paramsFitted_Individual.
        slopeStr = round(100*options.pfSlope);
        theData = load(fullfile(analysisDir,sprintf('%s_%s_ConstrainedSlope_%d.mat',options.subj,options.theIndicator,slopeStr)));
        thePFParams = theData.paramsFitted_Individual;
    end
else
    % In this case, parameters are in theData.paramsFittedAggregate.
    theData = load(fullfile(analysisDir,sprintf('%s_%s_%s.mat',options.subj,options.theIndicator,options.aggregateFitName)));
    thePFParams = theData.paramsFittedAggregate;
end

if isempty(theData)
    error('No fit data found');
end

%% Evaluate the PF
%
% First, round to the nearest .1 above the guess rate (code won't be happy if you
% try to evaluate the PF for prop seen below the lower asymptote)
%
% Note that the falsePosRate can vary with session if we don't correct for
% guessing, which might cause this to break if it's ever run over data
% aggregated over sessions.
PF = theData.PF;
propEvalStart = ceil(10*theData.sessionData{1}.falsePosRate)/10;
propSeen_Fit = propEvalStart:.1:.9;

modLevels_PF = nan(size(thePFParams,1), length(propSeen_Fit));
for n = 1:size(thePFParams,1)
    modLevels_PF(n,:) = 10.^PF(thePFParams(n,:), propSeen_Fit, 'inv');
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
        error('Check the logic in this code before relying on it.  Looks funny.');
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
theLim = options.theLim;
plot([-theLim theLim],[0 0],'k:','LineWidth',1);
plot([0 0],[-theLim theLim],'k:','LineWidth',1);
xlim([-theLim theLim]);
ylim([-theLim theLim]);
axis('square');
xlabel('Contrast 1')
ylabel('Contrast 2');
print(theEllipseFig,fullfile(analysisOutDir,sprintf('%s_%s_D%s_P%d_AllData.tiff', ...
    options.subj,options.computationalName,num2str(round(1000*options.defocusDiopters)),options.pupilDiam)), '-dtiff');

%% Fit ellipse
%
% The work of generating the points to plot is done by PointsOnEllipseQ.
% Can only do this if we have enough data.
nCirclePoints = 100;
circlePoints = UnitCircleGenerate(nCirclePoints);
if (size(theDataToFit, 2) > 4)
    [ellParams,A,Ainv,Q] = FitEllipseQ(theDataToFit,'lockAngleAt0',false,'errorScalar',errorScalar,'initialParams',[incrThresh(end) decrThresh(end) 45]);

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
    plot([-theLim theLim],[0 0],'k:','LineWidth',1);
    plot([0 0],[-theLim theLim],'k:','LineWidth',1);
    xlim([-theLim theLim]);
    ylim([-theLim theLim]);
    axis('square');
    xlabel('Contrast 1')
    ylabel('Contrast 2');
    print(theEllipseFig, fullfile(analysisOutDir,sprintf('%s_%s_D%s_P%d_EllipseAllData.tiff', ...
        options.subj,options.computationalName,num2str(round(1000*options.defocusDiopters)),options.pupilDiam)), '-dtiff');

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
    plot([-theLim theLim],[0 0],'k:','LineWidth',1);
    plot([0 0],[-theLim theLim],'k:','LineWidth',1);
    xlim([-theLim theLim]);
    ylim([-theLim theLim]);
    axis('square');
    xlabel('Contrast 1')
    ylabel('Contrast 2');
    print(theEllipseFig0, fullfile(analysisOutDir,sprintf('%s_%s_D%s_P%d_IncrAndDecrOnlyData.tiff', ...
        options.subj,options.computationalName,num2str(round(1000*options.defocusDiopters)),options.pupilDiam)), '-dtiff');
    plot(ellPoints0(1,:),ellPoints0(2,:),'r','LineWidth',2);
    print(theEllipseFig0, fullfile(analysisOutDir,sprintf('%s_%s_D%s_P%d_Ellipse0FitData.tiff', ...
        options.subj,options.computationalName,num2str(round(1000*options.defocusDiopters)),options.pupilDiam)), '-dtiff');

    % Add to plot of all data
    figure(theEllipseFig1);
    plot(ellPoints0(1,:),ellPoints0(2,:),'r','LineWidth',2);
    title(sprintf('Subject %s, Criterion %d%%',options.subj,round(100*propsSeen)));
    plot([-theLim theLim],[0 0],'k:','LineWidth',1);
    plot([0 0],[-theLim theLim],'k:','LineWidth',1);
    xlim([-theLim theLim]);
    ylim([-theLim theLim]);
    axis('square');
    xlabel('Contrast 1')
    ylabel('Contrast 2');
    print(theEllipseFig1, fullfile(analysisOutDir,sprintf('%s_%s_D%s_P%d_Ellipse0AllData.tiff', ...
        options.subj,options.computationalName,num2str(round(1000*options.defocusDiopters)),options.pupilDiam)), '-dtiff');
end

%% Get and fit ideal observer ellipse
%
% Directory stuff
baseProject = 'AOCompObserver';
compAnalysisBaseDir = getpref(baseProject,'analysisDir');
compAnalysisInDir = fullfile(compAnalysisBaseDir,sprintf('%s_%s_%d',options.computationalName,num2str(round(1000*options.defocusDiopters)),options.pupilDiam));
if (~exist(compAnalysisInDir ,'dir'))
    error('Computational observer not yet run for specified diopters of defocus');
end
compObserver = load(fullfile(compAnalysisInDir,sprintf('CompObserver_%s',options.computationalName)));
circlePointsFit(1,:) = cosd(stimAnglesFit);
circlePointsFit(2,:) = sind(stimAnglesFit);
compObserverEllData = PointsOnEllipseQ(compObserver.compFitQ,circlePointsFit);
for aa = 1:length(stimAnglesFit)
    dataRadii(aa) = norm(theDataToFit(:,aa));
    compRadii(aa) = norm(compObserverEllData(:,aa));
end
compFitFactor = compRadii'\dataRadii';
compObserverEll = PointsOnEllipseQ(compObserver.compFitQ,circlePoints)*compFitFactor;
theEllipseFig3 = figure; hold on
plot(compObserverEll(1,:),compObserverEll(2,:),'r','LineWidth',2);
title(sprintf('Subject %s, Criterion %d%%, Fit w/ Ideal, %0.2g D',options.subj,round(100*propsSeen),options.defocusDiopters));
fprintf('Ideal observer scale factor: %0.3g\n',compFitFactor);
plot([-theLim theLim],[0 0],'k:','LineWidth',1);
plot([0 0],[-theLim theLim],'k:','LineWidth',1);
xlim([-theLim theLim]);
ylim([-theLim theLim]);
axis('square');
xlabel('Contrast 1')
ylabel('Contrast 2');
print(theEllipseFig3, fullfile(analysisOutDir,sprintf('%s_%s_D%s_P%d_CompEllipse.tiff',...
    options.subj,options.computationalName,num2str(round(1000*options.defocusDiopters)),options.pupilDiam)), '-dtiff');
plot(theDataToFit(1,:),theDataToFit(2,:),[theColors(1) 'o'],'MarkerFaceColor',theColors(1),'MarkerSize',12);
print(theEllipseFig3, fullfile(analysisOutDir,sprintf('%s_%s_D%s_P%d_CompEllipseWithData.tiff', ...
    options.subj,options.computationalName,num2str(round(1000*options.defocusDiopters)),options.pupilDiam)), '-dtiff');

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
print(theIncrFig, fullfile(analysisOutDir,sprintf('%s_%s_D%s_P%d_SingleIncr.tiff', ...
    options.subj,options.computationalName,num2str(round(1000*options.defocusDiopters)),options.pupilDiam)), '-dtiff');

%% Both incremental stimuli
index = find(stimAnglesFit == 90);
if (iscell(theDataToFit))
    for pp = 1:length(propsSeen)
        plot(theDataToFit{pp}(1,index),theDataToFit{pp}(2,index),[theColors(pp) 'o'],'MarkerFaceColor',theColors(pp),'MarkerSize',12);
    end
else
    plot(theDataToFit(1,index),theDataToFit(2,index),[theColors(1) 'o'],'MarkerFaceColor',theColors(1),'MarkerSize',12);
end
print(theIncrFig, fullfile(analysisOutDir,sprintf('%s_%s_D%s_P%d_DoubleIncr.tiff', ...
    options.subj,options.computationalName,num2str(round(1000*options.defocusDiopters)),options.pupilDiam)), '-dtiff');
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
print(theIncrFig, fullfile(analysisOutDir,sprintf('%s_%s_D%s_P%d_DoubleIncrWithLine.tiff', ...
    options.subj,options.computationalName,num2str(round(1000*options.defocusDiopters)),options.pupilDiam)), '-dtiff');

%% Add incr-incr datum
index = find(stimAnglesFit == 45);
plot(theDataToFit(1,index),theDataToFit(2,index),[theColors(1) 'o'],'MarkerFaceColor',theColors(1),'MarkerSize',12);
print(theIncrFig, fullfile(analysisOutDir,sprintf('%s_%s_D%s_P%d_IncrWithDoubleIncrLine.tiff', ...
    options.subj,options.computationalName,num2str(round(1000*options.defocusDiopters)),options.pupilDiam)), '-dtiff');

%% Add incr-incr
figure(theIncrFig1);
index = find(stimAnglesFit == 45);
plot(theDataToFit(1,index),theDataToFit(2,index),[theColors(1) 'o'],'MarkerFaceColor',theColors(1),'MarkerSize',12);
print(theIncrFig1, fullfile(analysisOutDir,sprintf('%s_%s_D%s_P%d_Incr.tiff', ...
    options.subj,options.computationalName,num2str(round(1000*options.defocusDiopters)),options.pupilDiam)), '-dtiff');

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
print(theIncrFig1, fullfile(analysisOutDir,sprintf('%s_%s_D%s_P%d_IncrWithBestIncrLine.tiff', ...
    options.subj,options.computationalName,num2str(round(1000*options.defocusDiopters)),options.pupilDiam)), '-dtiff');

%% Add decr-decr points
neededAngles = [180 225 270];
index = [];
for aa = 1:length(neededAngles)
    index = [index find(stimAnglesFit == neededAngles(aa))];
end
plot(theDataToFit(1,index),theDataToFit(2,index),[theColors(1) 'o'],'MarkerFaceColor',theColors(1),'MarkerSize',12);
print(theIncrFig1, fullfile(analysisOutDir,sprintf('%s_%s_D%s_P%d_IncrAndDecrWithBestIncrLine.tiff', ...
    options.subj,options.computationalName,num2str(round(1000*options.defocusDiopters)),options.pupilDiam)), '-dtiff');

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
print(theIncrFig1, fullfile(analysisOutDir,sprintf('%s_%s_D%s_P%d_IncrAndDecrWithBestIncrAndDecrLines.tiff', ...
    options.subj,options.computationalName,num2str(round(1000*options.defocusDiopters)),options.pupilDiam)), '-dtiff');

%% Add the critical incr-decr points
index = find((stimAnglesFit > 90 & stimAnglesFit < 180) | (stimAnglesFit > 270 & stimAnglesFit < 360));
plot(theDataToFit(1,index),theDataToFit(2,index),[theColors(1) 'o'],'MarkerFaceColor',theColors(1),'MarkerSize',12);
print(theIncrFig1, fullfile(analysisOutDir,sprintf('%s_%s_D%s_P%d_AllDataWithBestIncrAndDecrLines.tiff', ...
    options.subj,options.computationalName,num2str(round(1000*options.defocusDiopters)),options.pupilDiam)), '-dtiff');

%% Save the data
close all;
save(fullfile(analysisOutDir,sprintf('%s_%s_D%s_P%d_ContourAnalysis', ...
    options.subj,options.computationalName,num2str(round(1000*options.defocusDiopters)),options.pupilDiam)));

