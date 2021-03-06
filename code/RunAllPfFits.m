% RunAllPFFits
%
% Run all the PF fits

fitPfsToIncrDecrData('subj','11043','dataDate','20200131','norm',false','corrGuess',true,'refl',false);
fitPfsToIncrDecrData('subj','11046','dataDate','20200131','norm',false','corrGuess',true,'refl',false);

fitPfsToIncrDecrData('subj','11043','dataDate','20200131','norm',false','corrGuess',true,'refl',true);
fitPfsToIncrDecrData('subj','11046','dataDate','20200131','norm',false','corrGuess',true,'refl',true);