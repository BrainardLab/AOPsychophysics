% RunAllPFFits
%
% Run all the PF fits

fitPfsToIncrDecrData('subj','11043','dataDate','20200131','subProject','IncrDecr1','condition','Separation_1', ...
    'norm',false','corrGuess',true,'refl',false);
fitPfsToIncrDecrData('subj','11046','dataDate','20200131','subProject','IncrDecr1','condition','Separation_1', ...
    'norm',false','corrGuess',true,'refl',false);

fitPfsToIncrDecrData('subj','11043','dataDate','20200131','subProject','IncrDecr1','condition','Separation_1', ...
    'norm',false','corrGuess',true,'refl',true);
fitPfsToIncrDecrData('subj','11046','dataDate','20200131','subProject','IncrDecr1','condition','Separation_1', ...
    'norm',false','corrGuess',true,'refl',true);


fitPfsToIncrDecrData('subj','11046','dataDate','20210914','subProject','IncrDecr2','condition','Size_1', ...
    'norm',false','corrGuess',true,'refl',false);
fitPfsToIncrDecrData('subj','11046','dataDate','20210914','subProject','IncrDecr2','condition','Size_1', ...
    'norm',false','corrGuess',true,'refl',true);

fitPfsToIncrDecrData('subj','11046','dataDate','20210914','subProject','IncrDecr2','condition','Size_2', ...
    'norm',false','corrGuess',true,'refl',false);
