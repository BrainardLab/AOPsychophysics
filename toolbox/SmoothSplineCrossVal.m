function [fitObj,smoothingParam] = SmoothSplineCrossVal(xVals,yVals,options)

arguments
    xVals
    yVals
    options.smoothingParamLow (1,1) = 0
    options.smoothingParamHigh (1,1) = 1
    options.nSmoothingParams (1,1) = 10
    options.nPartitions (1,1) = 10
    options.trainFraction (1,1) = 0.9
    options.plot (1,1) = true
end

% Get unique x values and sort
uniqueX = unique(xVals);
xUse = sort(uniqueX(:));

% Set up cross validation partitions
nPartitions = options.nPartitions;
trainFraction = options.trainFraction;
trainN = floor(trainFraction*length(xVals));
for pp = 1:nPartitions
    trainRaw = Shuffle(1:length(xVals));
    trainIndex = trainRaw(1:trainN);
    testIndex = trainRaw(trainN+1:end);
   
    xUseTrain{pp} = xVals(trainIndex);
    yUseTrain{pp} = yVals(trainIndex);
    xUseTest{pp} = xVals(testIndex);
    yUseTest{pp} = yVals(testIndex);
end

% Try out a bunch of smoothing parameters and cross validate
minSmooth = options.smoothingParamLow; 
maxSmooth = options.smoothingParamHigh;
nSmooth = options.nSmoothingParams;
smoothingParams = linspace(minSmooth,maxSmooth,nSmooth);
for ss = 1:nSmooth
    errSmooth(ss) = 0;
    for pp = 1:nPartitions
        smoothingParam = smoothingParams(ss);
        fitObj = fit(xUseTrain{pp},yUseTrain{pp},'smoothingspline','smoothingparam',smoothingParam);
        yFit = feval(fitObj,xUseTest{pp});
        errSmooth(ss) = errSmooth(ss) + sum((yUseTest{pp}-yFit).^2);
    end
end

% Set smoothing parameter and fit
[~,whichSmooth] = min(errSmooth);
smoothingParam = smoothingParams(whichSmooth);
fitObj = fit(xVals,yVals,'smoothingspline','smoothingparam',smoothingParam);

if (options.plot)
    figure; clf; hold on
    plot(smoothingParams,errSmooth,'ro','MarkerSize',12);
end

if (options.plot)
    figure; clf; hold on
    xSmooth = linspace(min(xUse),max(xUse),100)';
    ySmooth = feval(fitObj,xSmooth);
    plot(xVals,yVals,'ro','MarkerFaceColor','r','MarkerSize',12);
    plot(xSmooth,ySmooth,'r','LineWidth',2);
end