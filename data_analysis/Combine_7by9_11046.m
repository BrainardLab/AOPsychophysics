% Combine_7by9__11046
%
% High level script to combine 7x9 data across sessions
%
% This version for the 7x9 pixel spots, for
% subject 11046.

%% Clear
clear; close all;

%% Parameters
%
% Where's the analyzed data and precomputed computational observer output?
theSubject = '11046';
psychoProject = 'AOPsychophysics';
psychoBaseDir = getpref(psychoProject,'analysisDir');
dataFilename = '11046_incDecFits_Aggregated.mat';

% Computational observer parameters
compProject = 'AOCompObserver';
compBaseDir = getpref(compProject,'analysisDir');
computationalName = '7_9_0';
defocusDiopters = 0.10;
pupilDiam = 7;
compFilename = sprintf('%s_%s_D%s_P%d_ContourAnalysis.mat', ...
    theSubject,computationalName,num2str(round(1000*defocusDiopters)),pupilDiam);

% Normalization
SESSION_NORMALIZE = true;
INCRDECR_NORMALIZE = true;
ANGLE_AVERAGE = true;

% Specify sessions.
sessionNames = {'IncrDecr1', 'IncrDecr2', 'IncrDecr4', ...
    'IncrDecr5', 'IncrDecr5', 'IncrDecr5', 'IncrDecr5', 'IncrDecr5', ...
    'IncrDecr6', 'IncrDecr6', 'IncrDecr6', 'IncrDecr6', 'IncrDecr6', ...
    'IncrDecr6', 'IncrDecr6', 'IncrDecr6', 'IncrDecr6', 'IncrDecr6', ...
    'IncrDecr6', 'IncrDecr6', 'IncrDecr6', 'IncrDecr6', 'IncrDecr6', ...
    'IncrDecr6', 'IncrDecr6'};

idx = 1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20200131','Separation_1','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20210914','Size_1','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211026','Size_1','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211123','Size1_Sep0','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211123','Size1_Sep2','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211123','Size1_Sep4','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211123','Size1_Sep6','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211123','Size1_Sep8','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211217','Dir0_Sep0','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211217','Dir45_Sep0','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211217','Dir45_Sep4','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211217','Dir45_Sep8','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211217','Dir45_Sep12','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211217','Dir45_Sep16','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211217','Dir225_Sep0','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211217','Dir225_Sep4','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211217','Dir225_Sep8','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211217','Dir225_Sep12','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211217','Dir225_Sep16','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211217','Dir270_Sep0','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211217','Dir315_Sep0','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211217','Dir315_Sep4','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211217','Dir315_Sep8','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211217','Dir315_Sep12','notnorm_corrguess_norefl',dataFilename); idx = idx+1;
theFiles{idx} = fullfile(psychoBaseDir,sessionNames{idx},theSubject,'20211217','Dir315_Sep16','notnorm_corrguess_norefl',dataFilename); idx = idx+1;

% Other params
titleStr = '11046, 7x9';
theLim = 2;

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
    theDataFit = [theDataFit [theData{ii}.thresholdXData ; theData{ii}.thresholdYData] ];
    stimAnglesFit = [stimAnglesFit theData{ii}.stimAngleList];
    thresholdContrasts = [thresholdContrasts theData{ii}.thresholdContrasts];
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
if (any( (round(abs(stimAnglesFit-stimAnglesFitCheck)) ~= 0) & (round(abs(stimAnglesFit-stimAnglesFitCheck)) ~= 360) ))
    error('Angle inconsistency');
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
        sessionNormable(uu) = true;
        incrDecrNormable = true;
    else
        fprintf('No nomralization data available for session %d\n',uniqueSessions(uu));
        sessionNormable(uu) = false;
        theMarkers(uu) = 'o';
    end
end
if (incrDecrNormable)
    uNormIncrData = uNormIncrData/nNormingSessions;
    uNormDecrData = uNormDecrData/nNormingSessions;
    uNormData = uNormData/nNormingSessions;
end

% Normalize when possible if specified
SESSION_NORMALIZE = true;
for uu = 1:length(uniqueSessions)
    sessionIndices = find( (sessionNumbersFit == uniqueSessions(uu)) );
    if (SESSION_NORMALIZE)
        if (sessionNormable(uu))
            theDataFit = theDataFit*uNormData/normData(uu);
        else
            % Adjust plot marker if not possible
            theMarkers(uu) = '+';
        end
    end
end

% Take mean over angles
uniqueAngles = unique(stimAnglesFit);
for aa = 1:length(uniqueAngles)
    index = find( (stimAnglesFit == uniqueAngles(aa)) & (separationsFit == 0));
    thresholdContrast = mean(thresholdContrasts(index));
    theDataFitAvg(1,index) = thresholdContrast*cosd(uniqueAngles(aa));
    theDataFitAvg(2,index) = thresholdContrast*sind(uniqueAngles(aa));
end

% Normalize incr/decr if desired
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

%% Plot the data for 0 separation
theColors = ['r' 'g' 'b' 'c' 'y'];
theColor = 1;
theDataFig = figure; clf; hold on;

% Each session gets its own color, up to length of colors list
for uu = 1:length(uniqueSessions)
    theSessionIndex = find( (sessionNumbersFit == uniqueSessions(uu)) & (separationsFit == 0));
    plot(theDataFit(1,theSessionIndex),theDataFit(2,theSessionIndex),[theColors(theColor) theMarkers(uu)],'MarkerFaceColor',theColors(theColor),'MarkerSize',8);
    plot(theDataFitAvg(1,theSessionIndex),theDataFitAvg(2,theSessionIndex),['k' theMarkers(uu)],'MarkerFaceColor','k','MarkerSize',16);


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

%% Fit ellipse
errorScalar = 1000;
[ellParams,A,Ainv,ellQ] = FitEllipseQ(theDataFit,'lockAngleAt0',false,'errorScalar',errorScalar,'initialParams',[mean(thresholdContrasts) mean(thresholdContrasts) 0]);

% Get points on ellipse
nCirclePoints = 100;
circlePoints = UnitCircleGenerate(nCirclePoints);
ellPoints = PointsOnEllipseQ(ellQ,circlePoints);
%plot(ellPoints(1,:),ellPoints(2,:),'k:','LineWidth',1);

%% Read in ideal observer threshold contour and add to plot, after scaling
compAnalysisInDir = fullfile(compBaseDir,sprintf('%s_%s_%d',computationalName,num2str(round(1000*defocusDiopters)),pupilDiam));
if (~exist(compAnalysisInDir ,'dir'))
    error('Computational observer not yet run for specified diopters of defocus');
end
compObserver = load(fullfile(compAnalysisInDir,sprintf('CompObserver_%s',computationalName)));

% Scale comp observer ellipse to data
decrScaleFactors = 1;
for dd = 1:length(decrScaleFactors)
    circlePointsFit(1,:) = cosd(stimAnglesFit);
    circlePointsFit(2,:) = sind(stimAnglesFit);
    compObserverEllData = PointsOnEllipseQ(compObserver.compFitQ,circlePointsFit);
    for jj = 1:size(compObserverEllData,2)
        if (compObserverEllData(1,jj) < 0)
            compObserverEllData(1,jj) = compObserverEllData(1,jj)/decrScaleFactors(dd);
        end
        if (compObserverEllData(2,jj) < 0)
            compObserverEllData(2,jj) = compObserverEllData(2,jj)/decrScaleFactors(dd);
        end
    end
    for aa = 1:length(stimAnglesFit)
        dataRadii(aa) = vecnorm(theDataFit(:,aa));
        compRadii(aa) = vecnorm(compObserverEllData(:,aa));
    end
    compFitFactor = compRadii'\dataRadii';
    compObserverEllDataFit = compObserverEllData*compFitFactor;
    fitError(dd) = sqrt(sum((theDataFit(:)-compObserverEllDataFit(:)).^2));
end
[~,decrScaleIndex] = min(fitError);
decrScaleFactor = decrScaleFactors(decrScaleIndex);

compObserverEll = PointsOnEllipseQ(compObserver.compFitQ,circlePoints)*compFitFactor;
% for jj = 1:size(compObserverEllData,2)
%     if (compObserverEll(1,jj) < 0)
%         compObserverEll(1,jj) = compObserverEll(1,jj)/decrScaleFactor;
%     end
%     if (compObserverEll(2,jj) < 0)
%         compObserverEll(2,jj) = compObserverEll(2,jj)/decrScaleFactor;
%     end
% end

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
