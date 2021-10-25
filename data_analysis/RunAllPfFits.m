% RunAllPFFits
%
% Run all the PF fits.
%
% The refl option folds the data symmetrically over
% the up-down and down-up configurations.  Sometimes
% this is a cleaner way to look at it, and sometimes not.

%% 2020-01-31
%
% 7x9 pixel rectangular stimulus, 11043 
fitPfsToIncrDecrData('subj','11043','dataDate','20200131','subProject','IncrDecr1','condition','Separation_1', ...
    'norm',false','corrGuess',true,'refl',false);
fitPfsToIncrDecrData('subj','11046','dataDate','20200131','subProject','IncrDecr1','condition','Separation_1', ...
    'norm',false','corrGuess',true,'refl',false);

% 7x9 pixel rectangular stimulus, 11046
fitPfsToIncrDecrData('subj','11043','dataDate','20200131','subProject','IncrDecr1','condition','Separation_1', ...
    'norm',false','corrGuess',true,'refl',true);
fitPfsToIncrDecrData('subj','11046','dataDate','20200131','subProject','IncrDecr1','condition','Separation_1', ...
    'norm',false','corrGuess',true,'refl',true);

%% 2021-09-14
%
% 7x9 pixel rectangular stimulus, 11046
% This is a replication of the 2020-01-31 data.
fitPfsToIncrDecrData('subj','11046','dataDate','20210914','subProject','IncrDecr2','condition','Size_1', ...
    'norm',false','corrGuess',true,'refl',false);
fitPfsToIncrDecrData('subj','11046','dataDate','20210914','subProject','IncrDecr2','condition','Size_1', ...
    'norm',false','corrGuess',true,'refl',true);

% 5x7 pixel rectangular stimulus, 11046
%
% This was just one run and the stimulus levels were
% duplicated, which means it is a mess.  Data not 
% usable but motivated us to do it more carefully later
% on.
fitPfsToIncrDecrData('subj','11046','dataDate','20210914','subProject','IncrDecr2','condition','Size_2', ...
    'norm',false','corrGuess',true,'refl',false);
fitPfsToIncrDecrData('subj','11046','dataDate','20210914','subProject','IncrDecr2','condition','Size_2', ...
    'norm',false','corrGuess',true,'refl',false);

%% 2021-10-18
%
% 5x7 pixel rectangular stimulus, 11046
fitPfsToIncrDecrData('subj','11046','dataDate','20211018','subProject','IncrDecr2','condition','Size_2', ...
    'norm',false','corrGuess',true,'refl',false);
fitPfsToIncrDecrData('subj','11046','dataDate','20211018','subProject','IncrDecr2','condition','Size_2', ...
    'norm',false','corrGuess',true,'refl',true);

% 4x8 pixel rectangular stimulus, 11046
% 
% Just one run.  Don't see the incr-decr, and 
% hardly see anything else.  Too small given
% gamut limitations.
fitPfsToIncrDecrData('subj','11046','dataDate','20211018','subProject','IncrDecr2','condition','Size_3', ...
    'norm',false','corrGuess',true,'refl',true);
fitPfsToIncrDecrData('subj','11046','dataDate','20211018','subProject','IncrDecr2','condition','Size_3', ...
    'norm',false','corrGuess',true,'refl',false);