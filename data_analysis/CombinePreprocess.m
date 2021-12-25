%% CombinePreprocess
%
% Common preprocessing for various combination scripts.
% The scripts then differ in what they plot, etc.

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
        theData{ii}.stimAnglesList = ReflectAnglesAround45(theData{ii}.stimAnglesList);
    end
    stimAnglesFit = [stimAnglesFit CanonicalAngles(theData{ii}.stimAngleList)];

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
uniqueSeparations = unique(separationsFit);
for aa = 1:length(uniqueAngles)
    for ss = 1:length(uniqueSeparations)
        index = find( (stimAnglesFit == uniqueAngles(aa)) & (separationsFit == uniqueSeparations(ss)) );
        if (~isempty(index))
            thresholdContrast = mean(thresholdContrasts(index));
            theDataFitAvg(1,aa,ss) = thresholdContrast.*cosd(uniqueAngles(aa));
            theDataFitAvg(2,aa,ss) = thresholdContrast.*sind(uniqueAngles(aa));
        else
            theDataFitAvg(1,aa,ss) = NaN;
            theDataFitAvg(2,aa,ss) = NaN;
        end
    end
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