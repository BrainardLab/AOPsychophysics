% fitEllipsesToIncrDecrData
%
% Script to combine data from inc/dec experiments and fit PFs to various
% axes in the 2D stimulus space
%
% See also: FitEllipseQ, PointsOnEllipseQ, EllipsoidMatricesGenerate.

%% Housekeeping
clear; close all
baseProject = 'AOPsychophysics';
subProject = 'IncrDecr1';

%% Parameters
%
% By setting lockAngle to true, you get ellipses with axes aligned with
% x and y axis.
lockAngle = false;
if (lockAngle)
    lockAngleStr = 'lockAngle';
else
    lockAngleStr = 'noLockAngle';
end

%% Scale decrement thresholds to be the same as increments?
scaleDecr = true;
if (scaleDecr)
    scaleDecrStr = 'scaleDecr';
else
    scaleDecrStr = 'noScaleDecr';
end

%% Normalize?
%
% To normalize, or not to normalize? Select false to work with the raw
% modulations, true to normalize.
norm = false;
if (norm)
    normStr = 'norm';
else
    normStr = 'notnorm';
end

%% Correct for guessing?
corrGuess = true;
if (corrGuess)
    corrGuessStr = 'corrguess';
else
    corrGuessStr = 'notcorrguess';
end

%% Reflect data so that stim1 and stim2 are treated as symmetric?
reflIn = false;
reflOut = false;
if (reflIn)
    reflInStr = 'refl';
else
    reflInStr = 'noRefl';
end
if (reflOut)
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
norm = false;
if (norm)
    error('Need to implement analysis of fits from normalized modulations');
end

% Analyze fits with slopes constrained?
constrainedSlopeFits = true;
if (~constrainedSlopeFits)
    error('Need to implement analysis of unconstrained slope fits');
end

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

%% Set up directories.
%
% Note that we're reading the output of the program
% that fits psychometric functions, so the data comes from 'analysisDir'.
analysisSubDir = sprintf('%s_%s_%s',normStr,corrGuessStr,reflInStr);
analysisBaseDir = getpref(baseProject,'analysisDir');
analysisDir = fullfile(analysisBaseDir,subProject,subj,dataDate,'Separation_1',analysisSubDir);

%% Read output of psychometric fitting
theData = load(fullfile(analysisDir,sprintf('%s_incDecFits_ConstrainedSlope.mat',subj)));
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

%% Get incremet and decrement thresholds
for pp = 1:length(theData.stimAngleList)       
    if (theData.stimAngleList(pp) == 0)
        incrThresh = abs(cosd(theData.stimAngleList(pp)).*modLevels_PF(pp,:));
    elseif (theData.stimAngleList(pp) == 270)
        decrThresh = abs(sind(theData.stimAngleList(pp)).*modLevels_PF(pp,:));
    end      
end

%% Convert to x and y coordinates
%
% As we do this, we also reflect the symmetric data points,
% which makes the plots look better and also lets us
% make some didactically useful plots.
index = 1;
for pp = 1:length(theData.stimAngleList)
    % Get data in x,y form
    xPlot_Fit(index,:) = cosd(theData.stimAngleList(pp)).*modLevels_PF(pp,:);
    yPlot_Fit(index,:) = sind(theData.stimAngleList(pp)).*modLevels_PF(pp,:);
    angles_Fit(index) = theData.stimAngleList(pp);
    index = index+1;
    
    % This does the reflection.  Double up points on diagonal so that each
    % datum counts same number of times in the fit.
    if (reflOut)
        if (theData.stimAngleList(pp) == 0)
            xPlot_Fit(index,:) = cosd(90).*modLevels_PF(pp,:);
            yPlot_Fit(index,:) = sind(90).*modLevels_PF(pp,:);
            angles_Fit(index) = 90;
            index = index+1;
        elseif (theData.stimAngleList(pp) == 270)
            xPlot_Fit(index,:) = cosd(180).*modLevels_PF(pp,:);
            yPlot_Fit(index,:) = sind(180).*modLevels_PF(pp,:);
            angles_Fit(index) = 180;  
            index = index+1;
        elseif (theData.stimAngleList(pp) == 45)
            xPlot_Fit(index,:) = cosd(45).*modLevels_PF(pp,:);
            yPlot_Fit(index,:) = sind(45).*modLevels_PF(pp,:);
            angles_Fit(index) = 45;
            index = index+1;
        elseif (theData.stimAngleList(pp) == 225)
            xPlot_Fit(index,:) = cosd(225).*modLevels_PF(pp,:);
            yPlot_Fit(index,:) = sind(225).*modLevels_PF(pp,:);
            angles_Fit(index) = 225;
            index = index+1;
        elseif (theData.stimAngleList(pp) < 0)
            xPlot_Fit(index,:) = cosd(90-theData.stimAngleList(pp)).*modLevels_PF(pp,:);
            yPlot_Fit(index,:) = sind(90-theData.stimAngleList(pp)).*modLevels_PF(pp,:);
            angles_Fit(index) = 90-theData.stimAngleList(pp);
            index = index+1;
        end
    end
end

%% Make increments and decrements same scale, if desired
if (scaleDecr)
    incrDecrScaleFactor = incrThresh./decrThresh;
    for pp = 1:size(xPlot_Fit,1)
        for jj = 1:size(xPlot_Fit,2)
            if (xPlot_Fit(pp,jj) < 0)
                xPlot_Fit(pp,jj) = incrDecrScaleFactor(jj)*xPlot_Fit(pp,jj);
            end
            if (yPlot_Fit(pp,jj) < 0)
                yPlot_Fit(pp,jj) = incrDecrScaleFactor(jj)*yPlot_Fit(pp,jj);
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
% propsSeen = [0.50 0.70 0.90];
propsSeen = 0.70;
if (length(propsSeen) > 1)
    for pp = 1:length(propsSeen)
        whichCol = find(propSeen_Fit == propsSeen(pp));
        xData = xPlot_Fit(:,whichCol)';
        yData = yPlot_Fit(:,whichCol)';
        theDataToFit{pp} = [xData ; yData];
    end
else
    whichCol = find(propSeen_Fit == propsSeen);
    xData = xPlot_Fit(:,whichCol)';
    yData = yPlot_Fit(:,whichCol)';
    theDataToFit = [xData ; yData];
end

%% Plot the points
theColors = ['r' 'k' 'b' 'b' 'y' 'c'];
if (length(propsSeen) > length(theColors))
    error('Need to specify more plot colors, fewer propsSeen, or fix code to handle');
end
theEllipseFig = figure; clf; hold on;
if (iscell(theDataToFit))
    for pp = 1:length(propsSeen)
        plot(theDataToFit{pp}(1,:),theDataToFit{pp}(2,:),[theColors(pp) 'o'],'MarkerFaceColor',theColors(pp),'MarkerSize',12);
    end
else
    plot(theDataToFit(1,:),theDataToFit(2,:),[theColors(1) 'o'],'MarkerFaceColor',theColors(1),'MarkerSize',12); 
end

%% Fit.
[ellParams,A,Ainv,Q] = FitEllipseQ(theDataToFit,'lockAngleAt0',lockAngle,'errorScalar',errorScalar,'initialParams',[incrThresh(end) decrThresh(end) 45]);

%% Print ellipse parameters
%
% Give full length of axes, not radius equivalent (which is why it is 2 in
% the numerator of the numbers below, not 1).  Printout is a little
% different depending on whether angle is locked or not.
%
% There is a vestigal printout call if you only fit one ellipse and don't
% put it in a cell array.
if (lockAngle)
    if (iscell(theDataToFit))
        for cc = 1:length(theDataToFit)
            fprintf('x-axis length: %0.2f, y-axis length: %0.2f, major axis angle %0.1f deg\n',2/ellParams{cc}(1),2/ellParams{cc}(2),ellParams{cc}(3));
        end
    else
        fprintf('x-axis length: %0.2f, y-axis length: %0.2f, major axis angle %0.1f deg\n',2/ellParams(1),2/ellParams(2),ellParams(3));
    end
else
    if (iscell(theDataToFit))
        for cc = 1:length(theDataToFit)  
            fprintf('Major axis length: %0.2f, minor axis length: %0.2f, major axis angle %0.1f deg\n',2/ellParams{cc}(1),2/ellParams{cc}(2),ellParams{cc}(3));
        end
    else
        fprintf('Major axis length: %0.2f, minor axis length: %0.2f, major axis angle %0.1f deg\n',2/ellParams(1),2/ellParams(2),ellParams(3));   
    end
end

%% Add ellipses to plot
%
% The work of generating the points to plot is done by PointsOnEllipseQ.
nCirclePoints = 100;
circlePoints = UnitCircleGenerate(nCirclePoints);
ellPoints = PointsOnEllipseQ(Q,circlePoints);
if (iscell(ellPoints))
    titleStr = sprintf('Subject %s, critera',subj);
    for cc = 1:length(ellPoints)
        plot(ellPoints{cc}(1,:),ellPoints{cc}(2,:),theColors(cc),'LineWidth',2);
        titleStr = [titleStr sprintf(' %d%%',round(100*propsSeen(cc)))];
    end
    title(titleStr);
else
    plot(ellPoints(1,:),ellPoints(2,:),'r','LineWidth',2);
    title(sprintf('Subject %s, criterion %d%%',subj,round(100*propsSeen)));
end
theLim = 2;
plot([-theLim theLim],[0 0],'k:','LineWidth',1);
plot([0 0],[-theLim theLim],'k:','LineWidth',1);
xlim([-theLim theLim]);
ylim([-theLim theLim]);
axis('square');
xlabel('Contrast 1')
ylabel('Contrast 2');
print(theEllipseFig, fullfile(analysisDir,sprintf('%s_%s_%s_%s_%sEllipse.tiff', subj,lockAngleStr,scaleDecrStr,reflOutStr)), '-dtiff');

%% Build up some explanatory plots
%
% Just incremental stimulus one
theIncrFig = figure; clf; hold on
index = find(angles_Fit == 0);
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
print(theIncrFig, fullfile(analysisDir,sprintf('%s_%s_%s_%s_%sSingleIncr.tiff', subj,lockAngleStr,scaleDecrStr,reflOutStr)), '-dtiff');

%% Both incremental stimuli
index = find(angles_Fit == 90);
if (iscell(theDataToFit))
    for pp = 1:length(propsSeen)
        plot(theDataToFit{pp}(1,index),theDataToFit{pp}(2,index),[theColors(pp) 'o'],'MarkerFaceColor',theColors(pp),'MarkerSize',12);
    end
else
    plot(theDataToFit(1,index),theDataToFit(2,index),[theColors(1) 'o'],'MarkerFaceColor',theColors(1),'MarkerSize',12); 
end
print(theIncrFig, fullfile(analysisDir,sprintf('%s_%s_%s_%s_DoubleIncr.tiff', subj,lockAngleStr,scaleDecrStr,reflOutStr)), '-dtiff');
theIncrFig1 = theIncrFig.copy;
figure(theIncrFig);

%% Fit linear model to single increments and plot
neededAngles = [0 90];
index = [];
for aa = 1:length(neededAngles)
    index = [index find(angles_Fit == neededAngles(aa))];
end
x = theDataToFit(1,index)';
y = theDataToFit(2,index)';
lineParams = [x ones(size(x))]\y;
linePlotX = linspace(-theLim,theLim,100)';
linePlotY = [linePlotX ones(size(linePlotX))]*lineParams;
plot(linePlotX,linePlotY,'r','LineWidth',3);
print(theIncrFig, fullfile(analysisDir,sprintf('%s_%s_%s_%s_DoubleIncrWithLine.tiff', subj,lockAngleStr,scaleDecrStr,reflOutStr)), '-dtiff');

%% Add incr-incr datum
if (iscell(theDataToFit))
    for pp = 1:length(propsSeen)
        plot(theDataToFit{pp}(1,index),theDataToFit{pp}(2,index),[theColors(pp) 'o'],'MarkerFaceColor',theColors(pp),'MarkerSize',12);
    end
else
    plot(theDataToFit(1,index),theDataToFit(2,index),[theColors(1) 'o'],'MarkerFaceColor',theColors(1),'MarkerSize',12); 
end
print(theIncrFig1, fullfile(analysisDir,sprintf('%s_%s_%s_%s_IncrWithDoubleIncrLine.tiff', subj,lockAngleStr,scaleDecrStr,reflOutStr)), '-dtiff');

%% Add incr-incr
figure(theIncrFig1);
index = find(angles_Fit == 45);
if (iscell(theDataToFit))
    for pp = 1:length(propsSeen)
        plot(theDataToFit{pp}(1,index),theDataToFit{pp}(2,index),[theColors(pp) 'o'],'MarkerFaceColor',theColors(pp),'MarkerSize',12);
    end
else
    plot(theDataToFit(1,index),theDataToFit(2,index),[theColors(1) 'o'],'MarkerFaceColor',theColors(1),'MarkerSize',12); 
end
print(theIncrFig1, fullfile(analysisDir,sprintf('%s_%s_%s_%s_Incr.tiff', subj,lockAngleStr,scaleDecrStr,reflOutStr)), '-dtiff');

%% Fit linear model to single increments and plot
neededAngles = [0 45 90];
index = [];
for aa = 1:length(neededAngles)
    index = [index find(angles_Fit == neededAngles(aa))];
end
x = theDataToFit(1,index)';
y = theDataToFit(2,index)';
lineParams = [x ones(size(x))]\y;
linePlotX = linspace(-theLim,theLim,100)';
linePlotY = [linePlotX ones(size(linePlotX))]*lineParams;
plot(linePlotX,linePlotY,'r','LineWidth',3);
print(theIncrFig1, fullfile(analysisDir,sprintf('%s_%s_%s_%s_IncrWithBestIncrLine.tiff', subj,lockAngleStr,scaleDecrStr,reflOutStr)), '-dtiff');

%% Add decr-decr points
neededAngles = [180 225 270];
index = [];
for aa = 1:length(neededAngles)
    index = [index find(angles_Fit == neededAngles(aa))];
end
if (iscell(theDataToFit))
    for pp = 1:length(propsSeen)
        plot(theDataToFit{pp}(1,index),theDataToFit{pp}(2,index),[theColors(pp) 'o'],'MarkerFaceColor',theColors(pp),'MarkerSize',12);
    end
else
    plot(theDataToFit(1,index),theDataToFit(2,index),[theColors(1) 'o'],'MarkerFaceColor',theColors(1),'MarkerSize',12); 
end
print(theIncrFig1, fullfile(analysisDir,sprintf('%s_%s_%s_%s_IncrAndDecrWithBestIncrLine.tiff', subj,lockAngleStr,scaleDecrStr,reflOutStr)), '-dtiff');

%% Add decr-decr line
neededAngles = [180 225 270];
index = [];
for aa = 1:length(neededAngles)
    index = [index find(angles_Fit == neededAngles(aa))];
end
x = theDataToFit(1,index)';
y = theDataToFit(2,index)';
lineParams = [x ones(size(x))]\y;
linePlotX = linspace(-theLim,theLim,100)';
linePlotY = [linePlotX ones(size(linePlotX))]*lineParams;
plot(linePlotX,linePlotY,'r','LineWidth',3);
print(theIncrFig1, fullfile(analysisDir,sprintf('%s_%s_%s_%s_IncrAndDecrWithBestIncrAndDecrLines.tiff', subj,lockAngleStr,scaleDecrStr,reflOutStr)), '-dtiff');

%% Add the critical incr-decr points
index = find((angles_Fit > 90 & angles_Fit < 180) | (angles_Fit > 270 & angles_Fit < 360));
if (iscell(theDataToFit))
    for pp = 1:length(propsSeen)
        plot(theDataToFit{pp}(1,index),theDataToFit{pp}(2,index),[theColors(pp) 'o'],'MarkerFaceColor',theColors(pp),'MarkerSize',12);
    end
else
    plot(theDataToFit(1,index),theDataToFit(2,index),[theColors(1) 'o'],'MarkerFaceColor',theColors(1),'MarkerSize',12); 
end
print(theIncrFig1, fullfile(analysisDir,sprintf('%s_%s_%s_%s_AllDataWithBestIncrAndDecrLines.tiff', subj,lockAngleStr,scaleDecrStr,reflOutStr)), '-dtiff');
