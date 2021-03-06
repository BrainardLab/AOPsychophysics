% RunAllPFFits
%
% Run all the PF fits

fitThresholdContourIncrDecrData('subj','11043','dataDate','20200131','subProject','IncrDecr1','condition','Separation_1', ...
    'norm',false','corrGuess',true,'reflIn',false,'reflOut',false,'scaleDecr',false,'constraintedSlopeFits',false,'lockAngle',false);
fitThresholdContourIncrDecrData('subj','11046','dataDate','20200131','subProject','IncrDecr1','condition','Separation_1', ...
    'norm',false','corrGuess',true,'reflIn',false,'reflOut',false,'scaleDecr',false,'constraintedSlopeFits',false,'lockAngle',false);
