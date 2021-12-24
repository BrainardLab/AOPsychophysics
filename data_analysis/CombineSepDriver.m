% CombineSepDriver
%
% Description:
%    Script called after all parameters set up for a given
%    observer/condition.

%% Do the preprocessing
CombinePreprocess;

%% Read in ideal observer threshold contour and add to plot, after scaling
% compAnalysisInDir = fullfile(compBaseDir,sprintf('%s_%s_%d',computationalName,num2str(round(1000*defocusDiopters)),pupilDiam));
% if (~exist(compAnsepAnalysisInDir ,'dir'))
%     error('Computational observer not yet run for specified diopters of defocus');
% end
% compObserver = load(fullfile(compAnalysisInDir,sprintf('CompObserver_%s',computationalName)));
% 
% % Scale comp observer ellipse to data
% nCirclePoints = 100;
% circlePoints = UnitCircleGenerate(nCirclePoints);
% circlePointsFit(1,:) = cosd(stimAnglesFit);
% circlePointsFit(2,:) = sind(stimAnglesFit);
% compObserverEllData = PointsOnEllipseQ(compObserver.compFitQ,circlePointsFit);
% for jj = 1:size(compObserverEllData,2)
%     if (compObserverEllData(1,jj) < 0)
%         compObserverEllData(1,jj) = compObserverEllData(1,jj);
%     end
%     if (compObserverEllData(2,jj) < 0)
%         compObserverEllData(2,jj) = compObserverEllData(2,jj);
%     end
% end
% angleFitIndex = 1;
% for aa = 1:length(stimAnglesFit)
%     if (FIT_FIRSTTHIRDONLY)
%         if ( (stimAnglesFit(aa) >= 0 & stimAnglesFit(aa) <= 90) | ...
%                 (stimAnglesFit(aa) >= 180 & stimAnglesFit(aa) <= 270) )
%             dataRadii(angleFitIndex) = vecnorm(theDataFit(:,aa));
%             compRadii(angleFitIndex) = vecnorm(compObserverEllData(:,aa));
%             angleFitIndex = angleFitIndex+1;
%         end
%     else
%         dataRadii(angleFitIndex) = vecnorm(theDataFit(:,aa));
%         compRadii(angleFitIndex) = vecnorm(compObserverEllData(:,aa));
%         angleFitIndex = angleFitIndex+1;
%     end
% end
% compFitFactor = compRadii'\dataRadii';
% compObserverEllDataFit = compObserverEllData*compFitFactor;
% fitError = sqrt(sum((theDataFit(:)-compObserverEllDataFit(:)).^2));
% compObserverEll = PointsOnEllipseQ(compObserver.compFitQ,circlePoints)*compFitFactor;
% fprintf('Ideal observer scale factor: %0.3g\n',compFitFactor);

%% Plot the data 
theColors = ['r' 'g' 'b' 'c' 'y'];
theColor = 1;
theDataFig = figure; clf; 
set(gcf, 'Color', 'w', 'Units', 'inches', 'Position', [1 1 18 6]);
for ss = 1:length(anglesToAnalyze)
    % Set up subplot
    subplot(1,length(anglesToAnalyze),ss); hold on;

    % Plot average data versus separation
    angleIndex = find(uniqueAngles == anglesToAnalyze(ss));
    theFitDataAvgTemp = squeeze(theDataFitAvg(:,angleIndex,:));
    theresholdContrastsTemp = vecnorm(theFitDataAvgTemp);
    plot(uniqueSeparations,theresholdContrastsTemp,'o','Color',[0.4 0.4 0.4],'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize',16);

    % Find all data for this angle, and plot versus separation
    for uu = 1:length(uniqueSessions)
        index = find( (sessionNumbersFit == uniqueSessions(uu)) & (stimAnglesFit == anglesToAnalyze(ss)));
        theThresholdsToPlot = vecnorm(theDataFit(:,index));
        separationsToPlot = separationsFit(index);
        plot(separationsToPlot,theThresholdsToPlot,[theColors(theColor) 'o'],'MarkerFaceColor',theColors(theColor),'MarkerSize',8);

        % Cycle through colors
        theColor = theColor + 1;
        if (theColor > length(theColors))
            theColor = 1;
        end
    end
    xlim([0 20]); xlabel('Separation (Pixels)');
    ylim([0 theLim]); ylabel('Threshold Contrast');
    titleStr = { LiteralUnderscore(sprintf('%s_%s_D%s_P%d', ...
        theSubject,computationalName,num2str(round(1000*defocusDiopters)),pupilDiam)) , ...
        ; sprintf('Angle: %d',anglesToAnalyze(ss)) };
    title(titleStr);
end

% Plot ellipse at bottom 
% plot(compObserverEll(1,:),compObserverEll(2,:),'r','LineWidth',2);
% 
% % Then data averaged over each angle
% plot(theDataFitAvg(1,:),theDataFitAvg(2,:),'o','Color',[0.4 0.4 0.4],'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize',16);
% 
% % Each session gets its own color, up to length of colors list
% for uu = 1:length(uniqueSessions)
%     theSessionIndex = find( (sessionNumbersFit == uniqueSessions(uu)) & (separationsFit == 0));
%     plot(theDataFit(1,theSessionIndex),theDataFit(2,theSessionIndex),[theColors(theColor) theMarkers(uu)],'MarkerFaceColor',theColors(theColor),'MarkerSize',theMarkerSizes(uu));
% 
%     % Cycle through colors
%     theColor = theColor + 1;
%     if (theColor > length(theColors))
%         theColor = 1;
%     end
% end
% titleStr = LiteralUnderscore(sprintf('%s_%s_D%s_P%d', ...
%     theSubject,computationalName,num2str(round(1000*defocusDiopters)),pupilDiam));
% title(titleStr);
% plot([-theLim theLim],[0 0],'k:','LineWidth',1);
% plot([0 0],[-theLim theLim],'k:','LineWidth',1);
% xlim([-theLim theLim]);
% ylim([-theLim theLim]);
% axis('square');
% xlabel('Contrast 1')
% ylabel('Contrast 2');

%% Save figure
outDirname = 'aaCombinedSep';
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


