% CombineEllDriver
%
% Description:
%    Script called after all parameters set up for a given
%    observer/condition.

%% Load in data from each session
theDataFit = [];
stimAnglesFit = [];
thresholdContrasts = [];
separationsFit = [];
sessionNumbersFit = [];
for ii = 1:length(theFiles)
    warnState = warning('off','MATLAB:load:cannotInstantiateLoadedVariable');
    theData{ii} = load(theFiles{ii});
    warning(warnState);

    % Reflect stimulus angles around 45 degree line if desired.
    if (REFLECT)
        for aa = 1:length(theData{ii}.stimAngleList)
            if (theData{ii}.stimAngleList(aa) > 45 && theData{ii}.stimAngleList(aa) < 225)
                delta = theData{ii}.stimAngleList(aa) - 45;
                theData{ii}.stimAngleList(aa) = 45 - delta;
            else
                theData{ii}.stimAngleList(aa) = theData{ii}.stimAngleList(aa);
            end
        end
    end
    stimAnglesFit = [stimAnglesFit theData{ii}.stimAngleList];
    stimAnglesFit = CanonicalAngles(stimAnglesFit);

    % Deal with rouding error on angles if we reflected
    if (REFLECT)
        stimAnglesFit = MatchEntriesToTolerance(stimAnglesFit,angleTolerance);
    end

    thresholdContrasts = [thresholdContrasts theData{ii}.thresholdContrasts];
    xData = theData{ii}.thresholdContrasts.*cosd(theData{ii}.stimAngleList);
    yData = theData{ii}.thresholdContrasts.*sind(theData{ii}.stimAngleList);
    if (~REFLECT)
        if (max(abs(xData-theData{ii}.thresholdXData)) > 1e-10)
            error('Internal inconsistency');
        end
        if (max(abs(yData-theData{ii}.thresholdYData)) > 1e-10)
            error('Internal inconsistency');
        end
    end
    theDataFit = [theDataFit [xData ; yData] ];
    separationsFit = [separationsFit theData{ii}.stimSeparationPixels];
    sessionNumbersFit = [sessionNumbersFit theData{ii}.sessionNumbers];
end

% Consistency checks
[stimAnglesFitCheckRadians,thresholdContrastsFitCheck] = cart2pol(theDataFit(1,:),theDataFit(2,:));
stimAnglesFitCheck = rad2deg(stimAnglesFitCheckRadians);
thresholdContrastsFitCheck2 = vecnorm(theDataFit);
if (max(abs(thresholdContrasts-thresholdContrastsFitCheck)) > 1e-10)
    error('Threshold contrast inconsistency');
end
if (max(abs(thresholdContrasts-thresholdContrastsFitCheck2)) > 1e-10)
    error('Threshold contrast inconsistency via vecnorm');
end
if (~REFLECT)
    if (any( (round(abs(stimAnglesFit-stimAnglesFitCheck)) ~= 0) & (round(abs(stimAnglesFit-stimAnglesFitCheck)) ~= 360) ))
        error('Angle inconsistency');
    end
end

%% Find unique sessions
uniqueSessions = unique(sessionNumbersFit);

%% Normalize by sessions?
%
% Do we have normalization data for this session? If so, get it. 
nNormingSessions = 0;
uNormIncrData = 0;
uNormDecrData = 0;
uNormData = 0;
incrDecrNormable = false;
for uu = 1:length(uniqueSessions)
    pureIncrIndices{uu} = find( (sessionNumbersFit == uniqueSessions(uu)) & (stimAnglesFit == 0 | stimAnglesFit == 90) & separationsFit == 0);
    pureDecrIndices{uu} = find( (sessionNumbersFit == uniqueSessions(uu)) & (stimAnglesFit == 180 | stimAnglesFit == 270) & separationsFit == 0);
    if ((length(pureIncrIndices{uu}) + length(pureDecrIndices{uu}))== 4)
        normIncrData(uu) = mean(thresholdContrasts(pureIncrIndices{uu}));
        normDecrData(uu) = mean(thresholdContrasts(pureDecrIndices{uu}));
        normData(uu) = (normIncrData(uu)+normDecrData(uu))/2;

        uNormIncrData = uNormIncrData + normIncrData(uu);
        uNormDecrData = uNormDecrData + normDecrData(uu);
        uNormData = uNormData + normData(uu);
        nNormingSessions = nNormingSessions+1;

        fprintf('Normalization data available for session %d\n',uniqueSessions(uu));
        theMarkers(uu) = 'o';
        theMarkerSizes(uu) = 8;
        sessionNormable(uu) = true;
        incrDecrNormable = true;
    else
        fprintf('No nomralization data available for session %d\n',uniqueSessions(uu));
        sessionNormable(uu) = false;
        theMarkers(uu) = 'o';
        theMarkerSizes(uu) = 8;
    end
end
if (incrDecrNormable)
    uNormIncrData = uNormIncrData/nNormingSessions;
    uNormDecrData = uNormDecrData/nNormingSessions;
    uNormData = uNormData/nNormingSessions;
end

% Normalize when possible if specified
for uu = 1:length(uniqueSessions)
    sessionIndices = find( (sessionNumbersFit == uniqueSessions(uu)) );
    if (SESSION_NORMALIZE)
        if (sessionNormable(uu))
            theDataFit = theDataFit*uNormData/normData(uu);
            thresholdContrasts = vecnorm(theDataFit);
        else
            % Adjust plot marker if not possible
            theMarkers(uu) = '*';
            theMarkerSizes(uu) = 12;
        end
    end
end

% Take mean over angles
uniqueAngles = unique(stimAnglesFit);
for aa = 1:length(uniqueAngles)
    index = find( (stimAnglesFit == uniqueAngles(aa)) & (separationsFit == 0) );
    thresholdContrast = mean(thresholdContrasts(index));
    theDataFitAvg(1,aa) = thresholdContrast.*cosd(uniqueAngles(aa));
    theDataFitAvg(2,aa) = thresholdContrast.*sind(uniqueAngles(aa));
end

% Normalize incr/decr if desired
oldUniqueAngles = unique(CanonicalAngles(atan2d(theDataFit(2,:),theDataFit(1,:))));
if (INCRDECR_NORMALIZE)
    if (incrDecrNormable)
        index = find(theDataFit(1,:) < 0);
        theDataFit(1,index) = theDataFit(1,index)*uNormIncrData/uNormDecrData;%
        index = find(theDataFit(2,:) < 0);
        theDataFit(2,index) = theDataFit(2,index)*uNormIncrData/uNormDecrData;

        index = find(theDataFitAvg(1,:) < 0);
        theDataFitAvg(1,index) = theDataFitAvg(1,index)*uNormIncrData/uNormDecrData;
        index = find(theDataFitAvg(2,:) < 0);
        theDataFitAvg(2,index) = theDataFitAvg(2,index)*uNormIncrData/uNormDecrData;
    end
end
newUniqueAngles = unique(CanonicalAngles(atan2d(theDataFit(2,:),theDataFit(1,:))));

%% Read in ideal observer threshold contour and add to plot, after scaling
compAnalysisInDir = fullfile(compBaseDir,sprintf('%s_%s_%d',computationalName,num2str(round(1000*defocusDiopters)),pupilDiam));
if (~exist(compAnalysisInDir ,'dir'))
    error('Computational observer not yet run for specified diopters of defocus');
end
compObserver = load(fullfile(compAnalysisInDir,sprintf('CompObserver_%s',computationalName)));

% Scale comp observer ellipse to data
nCirclePoints = 100;
circlePoints = UnitCircleGenerate(nCirclePoints);
circlePointsFit(1,:) = cosd(stimAnglesFit);
circlePointsFit(2,:) = sind(stimAnglesFit);
compObserverEllData = PointsOnEllipseQ(compObserver.compFitQ,circlePointsFit);
for jj = 1:size(compObserverEllData,2)
    if (compObserverEllData(1,jj) < 0)
        compObserverEllData(1,jj) = compObserverEllData(1,jj);
    end
    if (compObserverEllData(2,jj) < 0)
        compObserverEllData(2,jj) = compObserverEllData(2,jj);
    end
end
angleFitIndex = 1;
for aa = 1:length(stimAnglesFit)
    if (FIT_FIRSTTHIRDONLY)
        if ( (stimAnglesFit(aa) >= 0 & stimAnglesFit(aa) <= 90) | ...
                (stimAnglesFit(aa) >= 180 & stimAnglesFit(aa) <= 270) )
            dataRadii(angleFitIndex) = vecnorm(theDataFit(:,aa));
            compRadii(angleFitIndex) = vecnorm(compObserverEllData(:,aa));
            angleFitIndex = angleFitIndex+1;
        end
    else
        dataRadii(angleFitIndex) = vecnorm(theDataFit(:,aa));
        compRadii(angleFitIndex) = vecnorm(compObserverEllData(:,aa));
        angleFitIndex = angleFitIndex+1;
    end
end
compFitFactor = compRadii'\dataRadii';
compObserverEllDataFit = compObserverEllData*compFitFactor;
fitError = sqrt(sum((theDataFit(:)-compObserverEllDataFit(:)).^2));
compObserverEll = PointsOnEllipseQ(compObserver.compFitQ,circlePoints)*compFitFactor;
fprintf('Ideal observer scale factor: %0.3g\n',compFitFactor);

%% Plot the data 
theColors = ['r' 'g' 'b' 'c' 'y'];
theColor = 1;
theDataFig = figure; clf; hold on;

% Plot ellipse at bottom 
plot(compObserverEll(1,:),compObserverEll(2,:),'r','LineWidth',2);

% Then data averaged over each angle
plot(theDataFitAvg(1,:),theDataFitAvg(2,:),'o','Color',[0.4 0.4 0.4],'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize',16);

% Each session gets its own color, up to length of colors list
for uu = 1:length(uniqueSessions)
    theSessionIndex = find( (sessionNumbersFit == uniqueSessions(uu)) & (separationsFit == 0));
    plot(theDataFit(1,theSessionIndex),theDataFit(2,theSessionIndex),[theColors(theColor) theMarkers(uu)],'MarkerFaceColor',theColors(theColor),'MarkerSize',theMarkerSizes(uu));

    % Cycle through colors
    theColor = theColor + 1;
    if (theColor > length(theColors))
        theColor = 1;
    end
end
titleStr = LiteralUnderscore(sprintf('%s_%s_D%s_P%d', ...
    theSubject,computationalName,num2str(round(1000*defocusDiopters)),pupilDiam));
title(titleStr);
plot([-theLim theLim],[0 0],'k:','LineWidth',1);
plot([0 0],[-theLim theLim],'k:','LineWidth',1);
xlim([-theLim theLim]);
ylim([-theLim theLim]);
axis('square');
xlabel('Contrast 1')
ylabel('Contrast 2');

%% Save figure
outDirname = 'aaCombinedEll';
outputPath = fullfile(psychoBaseDir,outDirname);
if (~exist(outputPath,'dir'))
    mkdir(outputPath);
end
outputFilename = sprintf('%s_%s_D%s_P%d_Ellipse.tiff', ...
    theSubject,computationalName,num2str(round(1000*defocusDiopters)),pupilDiam);
print(theDataFig, fullfile(outputPath,outputFilename), '-dtiff');

%% Fit ellipse and add to plot
% errorScalar = 1000;
% [ellParams,A,Ainv,ellQ] = FitEllipseQ(theDataFit,'lockAngleAt0',false,'errorScalar',errorScalar,'initialParams',[mean(thresholdContrasts) mean(thresholdContrasts) 0]);
% nCirclePoints = 100;
% circlePoints = UnitCircleGenerate(nCirclePoints);
% ellPoints = PointsOnEllipseQ(ellQ,circlePoints);
% plot(ellPoints(1,:),ellPoints(2,:),'k:','LineWidth',1);


