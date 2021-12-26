% CombineSepDriver
%
% Description:
%    Script called after all parameters set up for a given
%    observer/condition.
%
%    To focus on effect of separation, this combines across up-down/down-up
%    pairs.

%% Do the preprocessing
CombinePreprocess;

%% Plot the data 
theColors = ['r' 'g' 'b' 'c' 'k' 'y'];
theDataFig = figure; clf; 
set(gcf, 'Color', 'w', 'Units', 'inches', 'Position', [1 1 18 6]);
for aa = 1:length(anglesToAnalyze)
    % Set up subplot
    subplot(1,length(anglesToAnalyze),aa); hold on;
    theColor = 1;

    % Plot average data versus separation.  Because we're combining over
    % symmetric angles, need to take the average of the contrasts over the
    % second dimension here.
    angleIndex = find(ReflectAnglesAround45(uniqueAngles) == anglesToAnalyze(aa));
    if (length(angleIndex) == 1)
        theFitDataAvgTemp = squeeze(theDataFitAvg(:,angleIndex,:));
        thresholdContrastsAvg{aa} = vecnorm(theFitDataAvgTemp);
    else
        theFitDataAvgTemp = squeeze(theDataFitAvg(:,angleIndex,:));
        thresholdContrastsAvg{aa} = squeeze(mean(vecnorm(theFitDataAvgTemp),2));
    end
    plot(uniqueSeparations*minPerPixel,thresholdContrastsAvg{aa},'o','Color',[0.4 0.4 0.4],'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize',16);

    % Find all data for this angle, and plot versus separation
    for uu = 1:length(uniqueSessions)
        index = find( (sessionNumbersFit == uniqueSessions(uu)) & (ReflectAnglesAround45(stimAnglesFit) == anglesToAnalyze(aa)));
        theThresholdsToPlot = vecnorm(theDataFit(:,index));
        separationsToPlot = separationsFit(index);
        plot(separationsToPlot*minPerPixel,theThresholdsToPlot,[theColors(theColor) 'o'],'MarkerFaceColor',theColors(theColor),'MarkerSize',6);

        % Cycle through colors
        theColor = theColor + 1;
        if (theColor > length(theColors))
            theColor = 1;
        end
    end
    xlim([0 20*minPerPixel]); xlabel('Separation (arcmin)');
    ylim([0 theLim]); ylabel('Threshold Contrast');
    titleStr = { LiteralUnderscore(sprintf('%s_%s_D%s_P%d', ...
        theSubject,computationalName,num2str(round(1000*defocusDiopters)),pupilDiam)) , ...
        ; sprintf('Angle: %d',anglesToAnalyze(aa)) };
    title(titleStr);
end

%% Read in ideal observer threshold contour and scale to data
for ss = 1:length(compSeparations)
    compAnalysisInDir = fullfile(compBaseDir,sprintf('%s_%d_%s_%d',computationalName,compSeparations(ss),num2str(round(1000*defocusDiopters)),pupilDiam));
    if (~exist(compAnalysisInDir ,'dir'))
        error('Computational observer not yet run for specified diopters of defocus');
    end
    compObserver{ss} = load(fullfile(compAnalysisInDir,sprintf('CompObserver_%s_%d',computationalName,compSeparations(ss))));
end

% For each angle to analyze, find comp observer threshold for each
% separation and set up regression to average data.
regressData = [];
regressComp = [];
for aa = 1:length(anglesToAnalyze)

    % Optional only fit to first and third quadrants.  This logic skips
    % angles in second and fourth quadrant if the flag is true.
    regressIt = false;
    if (FIT_FIRSTTHIRDONLY)
        if ( (anglesToAnalyze(aa) >= 0 & anglesToAnalyze(aa) <= 90) | ...
                (anglesToAnalyze(aa) >= 180 & anglesToAnalyze(aa) <= 270) )
            regressIt = true;
        end
    else
        regressIt = true;
    end

    % Get comp observer points
    circlePointsFit(1,:) = cosd(anglesToAnalyze(aa));
    circlePointsFit(2,:) = sind(anglesToAnalyze(aa));
    index = [];
    for ss = 1:length(compSeparations)
        tempData = PointsOnEllipseQ(compObserver{ss}.compFitQ,circlePointsFit);
        compObserverSepData(aa,ss) = vecnorm(tempData);
        idx = find(compSeparations(ss) == uniqueSeparations);
        index = [index idx];
    end

    % Add to regression data, optionally
    if (regressIt)
        regressData = [regressData ; thresholdContrastsAvg{aa}'];
        regressComp = [regressComp ; compObserverSepData(aa,index)'];
    end
end

% Scale comp observer functions to observer data
compFitFactor = regressComp(~isnan(regressData))\regressData(~isnan(regressData));
fprintf('Ideal observer scale factor: %0.3g\n',compFitFactor);
compObserverSepDataScaled = compFitFactor*compObserverSepData;

% Add scaled comp observer to plot
sepSmoothPoints = linspace(min(compSeparations),max(compSeparations),100);
for aa = 1:length(anglesToAnalyze)
    % Smooth spline through the computational points
    smoothingParameter = 0.96;
    fSpline = fit(compSeparations'*minPerPixel,compObserverSepDataScaled(aa,:)','smoothingspline','SmoothingParam',smoothingParameter);
    splinePredictions = feval(fSpline,sepSmoothPoints);

    % Set up subplot
    subplot(1,length(anglesToAnalyze),aa); hold on;
    plot(compSeparations*minPerPixel,compObserverSepDataScaled(aa,:),'r*','MarkerSize',6);
    plot(sepSmoothPoints,splinePredictions,'r','LineWidth',2);
end

%% Save figure
outDirname = 'aaCombinedSep';
outputPath = fullfile(psychoBaseDir,outDirname);
if (~exist(outputPath,'dir'))
    mkdir(outputPath);
end
outputFilename = sprintf('%s_%s_D%s_P%d_Ellipse.tiff', ...
    theSubject,computationalName,num2str(round(1000*defocusDiopters)),pupilDiam);
print(theDataFig, fullfile(outputPath,outputFilename), '-dtiff');
