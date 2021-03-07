% RunAllPFFits
%
% Run all the PF fits

% Basic fits, no decr scaling
fitThresholdContourIncrDecrData('subj','11043','dataDate','20200131','subProject','IncrDecr1','condition','Separation_1', ...
    'norm',false','corrGuess',true,'reflIn',false,'reflOut',false,'scaleDecr',false,'constrainedSlopeFits',true,'lockAngle',false,'lockSlope',true);
fitThresholdContourIncrDecrData('subj','11046','dataDate','20200131','subProject','IncrDecr1','condition','Separation_1', ...
    'norm',false','corrGuess',true,'reflIn',false,'reflOut',false,'scaleDecr',false,'constrainedSlopeFits',true,'lockAngle',false,'lockSlope',true);

% Basic fits, with decr scaling
fitThresholdContourIncrDecrData('subj','11043','dataDate','20200131','subProject','IncrDecr1','condition','Separation_1', ...
    'norm',false','corrGuess',true,'reflIn',false,'reflOut',false,'scaleDecr',true,'constrainedSlopeFits',true,'lockAngle',false,'lockSlope',true);
fitThresholdContourIncrDecrData('subj','11046','dataDate','20200131','subProject','IncrDecr1','condition','Separation_1', ...
    'norm',false','corrGuess',true,'reflIn',false,'reflOut',false,'scaleDecr',true,'constrainedSlopeFits',true,'lockAngle',false,'lockSlope',true);

% Reflect out, with decr scaling
fitThresholdContourIncrDecrData('subj','11043','dataDate','20200131','subProject','IncrDecr1','condition','Separation_1', ...
    'norm',false','corrGuess',true,'reflIn',false,'reflOut',true,'scaleDecr',true,'constrainedSlopeFits',true,'lockAngle',false,'lockSlope',true);
fitThresholdContourIncrDecrData('subj','11046','dataDate','20200131','subProject','IncrDecr1','condition','Separation_1', ...
    'norm',false','corrGuess',true,'reflIn',false,'reflOut',true,'scaleDecr',true,'constrainedSlopeFits',true,'lockAngle',false,'lockSlope',true);


