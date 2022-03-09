function [fitObj,smoothingParam] = SmoothSplineCrossVal(xVals,yVals,options)

arguments
    xVals
    yVals
    options.smoothingParam (1,1) = 0.9
    options.plot (1,1) = true;
end

% Get unique x values and sort
uniqueX = unique(xVals);
xUse = sort(uniqueX(:));

% Set up cross validation partitions
nPartitions = 10;
trainFraction = 0.9;
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


% for uu = 1:length(xUse)
%     index = find(xVals == xUse(uu));
%     yUse(uu) = mean(yVals(index));
% end

minSmooth = 0; 
maxSmooth = 0.05;
nSmooth = 100;
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
smoothingParam = 0.001;
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