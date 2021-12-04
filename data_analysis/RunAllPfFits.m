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
clear; close all;
fitPfsToIncrDecrData('subj','11043','dataDate','20200131','subProject','IncrDecr1','condition','Separation_1', ...
    'norm',false','corrGuess',true,'refl',false);
fitPfsToIncrDecrData('subj','11046','dataDate','20200131','subProject','IncrDecr1','condition','Separation_1', ...
    'norm',false','corrGuess',true,'refl',false);

% 7x9 pixel rectangular stimulus, 11046
clear; close all;
fitPfsToIncrDecrData('subj','11043','dataDate','20200131','subProject','IncrDecr1','condition','Separation_1', ...
    'norm',false','corrGuess',true,'refl',true);
fitPfsToIncrDecrData('subj','11046','dataDate','20200131','subProject','IncrDecr1','condition','Separation_1', ...
    'norm',false','corrGuess',true,'refl',true);

%% 2021-09-14
%
% 7x9 pixel rectangular stimulus, 11046
% This is a replication of the 2020-01-31 data.
clear; close all;
fitPfsToIncrDecrData('subj','11046','dataDate','20210914','subProject','IncrDecr2','condition','Size_1', ...
    'norm',false','corrGuess',true,'refl',false);
fitPfsToIncrDecrData('subj','11046','dataDate','20210914','subProject','IncrDecr2','condition','Size_1', ...
    'norm',false','corrGuess',true,'refl',true);

% 5x7 pixel rectangular stimulus, 11046
%
% This was just one run and the stimulus levels were duplicated, which
% means it is a mess.  Data not usable but motivated us to do it more
% carefully later on.
%
% To run this, you have to uncomment some kluges in fitPfsToIncrDecrData
% which are currently commented out.  Those kluges massage the data a
% little so it looks like there was a fully crossed stimulus design. Search
% on string 9/14/21 to find the relevant code sections.
%
% fitPfsToIncrDecrData('subj','11046','dataDate','20210914','subProject','IncrDecr2','condition','Size_2', ...
%     'norm',false','corrGuess',true,'refl',false);
% fitPfsToIncrDecrData('subj','11046','dataDate','20210914','subProject','IncrDecr2','condition','Size_2', ...
%     'norm',false','corrGuess',true,'refl',false);

%% 2021-10-18
%
% 5x7 pixel rectangular stimulus, 11046
% clear; close all;
% fitPfsToIncrDecrData('subj','11046','dataDate','20211018','subProject','IncrDecr3','condition','Size_2', ...
%     'norm',false','corrGuess',true,'refl',false);
% fitPfsToIncrDecrData('subj','11046','dataDate','20211018','subProject','IncrDecr3','condition','Size_2', ...
%     'norm',false','corrGuess',true,'refl',true);

% 4x8 pixel rectangular stimulus, 11046
% 
% Just one run.  Don't see the incr-decr, and hardly see anything else.
% Too small given gamut limitations.
%
% For some reason I haven't tracked down, the 'refl',false case crashes,
% complaing about violation of the same number of levels in each direction.
% I have not tracked that down.
%
% fitPfsToIncrDecrData('subj','11046','dataDate','20211018','subProject','IncrDecr2','condition','Size_3', ...
%     'norm',false','corrGuess',true,'refl',true);
% fitPfsToIncrDecrData('subj','11046','dataDate','20211018','subProject','IncrDecr2','condition','Size_3', ...
%     'norm',false','corrGuess',true,'refl',false);

%% 2021-10-26
%
% 5x7 pixel rectangular stimulus, 11046. 4th quadrant conditions
% clear; close all;
% fitPfsToIncrDecrData('subj','11046','dataDate','20211026','subProject','IncrDecr4','condition','Size_2', ...
%     'norm',false','corrGuess',true,'refl',false);
% fitPfsToIncrDecrData('subj','11046','dataDate','20211026','subProject','IncrDecr4','condition','Size_2', ...
%     'norm',false','corrGuess',true,'refl',true);

% 6x8 pixel rectangular stimulus, 11046. 4th quadrant conditions
clear; close all;
fitPfsToIncrDecrData('subj','11046','dataDate','20211026','subProject','IncrDecr4','condition','Size_4', ...
    'norm',false','corrGuess',true,'refl',false);
fitPfsToIncrDecrData('subj','11046','dataDate','20211026','subProject','IncrDecr4','condition','Size_4', ...
    'norm',false','corrGuess',true,'refl',true);

% 7x9 pixel rectangular stimulus, 11046. 4th quadrant conditions
clear; close all;
fitPfsToIncrDecrData('subj','11046','dataDate','20211026','subProject','IncrDecr4','condition','Size_1', ...
    'norm',false','corrGuess',true,'refl',false);
fitPfsToIncrDecrData('subj','11046','dataDate','20211026','subProject','IncrDecr4','condition','Size_1', ...
    'norm',false','corrGuess',true,'refl',true);

%% 2021-11-23
%
% 7x9 pixel rectangular stimulus, 11046. 315 degrees, 0, 2, 4, 6 8 pixel
% separation
clear; close all;
fitPfsToIncrDecrData('subj','11046','dataDate','20211123','subProject','IncrDecr5','condition','Size1_Sep0', ...
    'norm',false','corrGuess',true,'refl',true);
fitPfsToIncrDecrData('subj','11046','dataDate','20211123','subProject','IncrDecr5','condition','Size1_Sep2', ...
    'norm',false','corrGuess',true,'refl',true);
fitPfsToIncrDecrData('subj','11046','dataDate','20211123','subProject','IncrDecr5','condition','Size1_Sep4', ...
    'norm',false','corrGuess',true,'refl',true);
fitPfsToIncrDecrData('subj','11046','dataDate','20211123','subProject','IncrDecr5','condition','Size1_Sep6', ...
    'norm',false','corrGuess',true,'refl',true);
fitPfsToIncrDecrData('subj','11046','dataDate','20211123','subProject','IncrDecr5','condition','Size1_Sep8', ...
    'norm',false','corrGuess',true,'refl',true);

% Fix slope to something reasonable
theSlope = 12;
fitPfsToIncrDecrData('subj','11046','dataDate','20211123','subProject','IncrDecr5','condition','Size1_Sep0', ...
    'norm',false','corrGuess',true,'refl',true,'fixedSlope',theSlope);
fitPfsToIncrDecrData('subj','11046','dataDate','20211123','subProject','IncrDecr5','condition','Size1_Sep0', ...
    'norm',false','corrGuess',true,'refl',false,'fixedSlope',theSlope);
fitPfsToIncrDecrData('subj','11046','dataDate','20211123','subProject','IncrDecr5','condition','Size1_Sep2', ...
    'norm',false','corrGuess',true,'refl',true,'fixedSlope',theSlope);
fitPfsToIncrDecrData('subj','11046','dataDate','20211123','subProject','IncrDecr5','condition','Size1_Sep4', ...
    'norm',false','corrGuess',true,'refl',true,'fixedSlope',theSlope);
fitPfsToIncrDecrData('subj','11046','dataDate','20211123','subProject','IncrDecr5','condition','Size1_Sep6', ...
    'norm',false','corrGuess',true,'refl',true,'fixedSlope',theSlope);
fitPfsToIncrDecrData('subj','11046','dataDate','20211123','subProject','IncrDecr5','condition','Size1_Sep8', ...
    'norm',false','corrGuess',true,'refl',true,'fixedSlope',theSlope);


% 6x8 pixel rectangular stimulus, 11046. 315 degrees, 0 and 4 pixel
% separation
clear; close all;
fitPfsToIncrDecrData('subj','11046','dataDate','20211123','subProject','IncrDecr5','condition','Size2_Sep0', ...
    'norm',false','corrGuess',true,'refl',true);
fitPfsToIncrDecrData('subj','11046','dataDate','20211123','subProject','IncrDecr5','condition','Size2_Sep4', ...
    'norm',false','corrGuess',true,'refl',true);

% Fix slope
theSlope = 12;
fitPfsToIncrDecrData('subj','11046','dataDate','20211123','subProject','IncrDecr5','condition','Size2_Sep0', ...
    'norm',false','corrGuess',true,'refl',true,'fixedSlope',theSlope);
fitPfsToIncrDecrData('subj','11046','dataDate','20211123','subProject','IncrDecr5','condition','Size2_Sep4', ...
    'norm',false','corrGuess',true,'refl',true,'fixedSlope',theSlope);

% Final close up.
clear; close all;