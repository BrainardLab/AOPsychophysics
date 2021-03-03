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
% By setting lockAngleAt0 to true, you get ellipses with axes aligned with
% x and y axis.
lockAngleAt0 = false;

%% Correct for guessing
correctForGuessing = true;

%% Reflect data to mirror symmetry of two dots in the experiment
reflect = true;

%% Scale increment and decrement thresholds to be the same
scaleIncrDecr = true;

%% If fit gets stuck, try futzing with the value of errorScalar.
errorScalar = 1000;

%% Which fits?
%
% To normalize, or not to normalize? Select false to work with the fits
% from the raw modulations, true for those from the normalized ones.
normFlag = false;
if (normFlag)
    error('Need to implement analysis of fits from normalized modulations');
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
analysisBaseDir = getpref(baseProject,'analysisDir');
analysisDir = fullfile(analysisBaseDir,subProject,subj,dataDate,'Separation_1');

%% Read output of psychometric fitting
theData = load(fullfile(analysisDir,sprintf('%s_%d_%d_incDecFits_ConstrainedSlope.mat',subj,normFlag,correctForGuessing)));
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
for ii = 1:length(theData.stimAngleList)       
    if (theData.stimAngleList(ii) == 0)
        incrThresh = abs(cosd(theData.stimAngleList(ii)).*modLevels_PF(ii,:));
    elseif (theData.stimAngleList(ii) == 270)
        decrThresh = abs(sind(theData.stimAngleList(ii)).*modLevels_PF(ii,:));
    end      
end

%% Convert to x and y coordinates
%
% As we do this, we also reflect the symmetric data points,
% which makes the plots look better and also lets us
% make some didactically useful plots.
index = 1;
for ii = 1:length(theData.stimAngleList)
    % Get data in x,y form
    xPlot_Fit(index,:) = cosd(theData.stimAngleList(ii)).*modLevels_PF(ii,:);
    yPlot_Fit(index,:) = sind(theData.stimAngleList(ii)).*modLevels_PF(ii,:);
    index = index+1;
    
    % This does the reflection.  Double up points on diagonal so that each
    % datum counts same number of times in the fit.
    if (reflect)
        if (theData.stimAngleList(ii) == 0)
            xPlot_Fit(index,:) = cosd(90).*modLevels_PF(ii,:);
            yPlot_Fit(index,:) = sind(90).*modLevels_PF(ii,:);
            index = index+1;
        elseif (theData.stimAngleList(ii) == 270)
            xPlot_Fit(index,:) = cosd(180).*modLevels_PF(ii,:);
            yPlot_Fit(index,:) = sind(180).*modLevels_PF(ii,:);
            index = index+1;
        elseif (theData.stimAngleList(ii) == 45)
            xPlot_Fit(index,:) = cosd(45).*modLevels_PF(ii,:);
            yPlot_Fit(index,:) = sind(45).*modLevels_PF(ii,:);
            index = index+1;
        elseif (theData.stimAngleList(ii) == 225)
            xPlot_Fit(index,:) = cosd(225).*modLevels_PF(ii,:);
            yPlot_Fit(index,:) = sind(225).*modLevels_PF(ii,:);
            index = index+1;
        elseif (theData.stimAngleList(ii) < 0)
            xPlot_Fit(index,:) = cosd(90-theData.stimAngleList(ii)).*modLevels_PF(ii,:);
            yPlot_Fit(index,:) = sind(90-theData.stimAngleList(ii)).*modLevels_PF(ii,:);
            index = index+1;
        end
    end
end

%% Make increments and decrements same scale, if desired
if (scaleIncrDecr)
    incrDecrScaleFactor = incrThresh./decrThresh;
    for ii = 1:size(xPlot_Fit,1)
        for jj = 1:size(xPlot_Fit,2)
            if (xPlot_Fit(ii,jj) < 0)
                xPlot_Fit(ii,jj) = incrDecrScaleFactor(jj)*xPlot_Fit(ii,jj);
            end
            if (yPlot_Fit(ii,jj) < 0)
                yPlot_Fit(ii,jj) = incrDecrScaleFactor(jj)*yPlot_Fit(ii,jj);
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
    for ii = 1:length(propsSeen)
        whichCol = find(propSeen_Fit == propsSeen(ii));
        xData = xPlot_Fit(:,whichCol)';
        yData = yPlot_Fit(:,whichCol)';
        theDataToFit{ii} = [xData ; yData];
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
theFig = figure; clf; hold on;
if (iscell(theDataToFit))
    for ii = 1:length(propsSeen)
        plot(theDataToFit{ii}(1,:),theDataToFit{ii}(2,:),[theColors(ii) 'o'],'MarkerFaceColor',theColors(ii),'MarkerSize',12);
    end
else
    plot(theDataToFit(1,:),theDataToFit(2,:),[theColors(1) 'o'],'MarkerFaceColor',theColors(1),'MarkerSize',12); 
end

theLim = 1;
plot([-theLim theLim],[0 0],'k:','LineWidth',1);
plot([0 0],[-theLim theLim],'k:','LineWidth',1);
xlim([-theLim theLim]);
ylim([-theLim theLim]);
axis('square');

%% Fit.
[ellParams,A,Ainv,Q] = FitEllipseQ(theDataToFit,'lockAngleAt0',lockAngleAt0,'errorScalar',errorScalar,'initialParams',[incrThresh(end) decrThresh(end) 45]);

%% Print ellipse parameters
%
% Give full length of axes, not radius equivalent (which is why it is 2 in
% the numerator of the numbers below, not 1).  Printout is a little
% different depending on whether angle is locked or not.
%
% There is a vestigal printout call if you only fit one ellipse and don't
% put it in a cell array.
if (lockAngleAt0)
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
xlabel('Stimulus 1')
ylabel('Stimulus 2');
print(gcf, fullfile(analysisDir,sprintf('%s_%d_%d_Ellipse.tiff', subj,normFlag,correctForGuessing)), '-dtiff');

%% Save fit data to mat file
% save(fullfile(analysisDir,sprintf('%s_%d_incDecFits_ConstrainedSlope.mat', subj,normFlag)), 'stimAngleList', 'falsePosProp', 'paramsFitted_Individual', 'paramsFitted_Multi', 'PF');
% print(gcf, fullfile(analysisDir,sprintf('%s_%d_incDecFits_ConstrainedSlope.png', subj,normFlag)), '-dpng2');
