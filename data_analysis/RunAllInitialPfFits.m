% RunAllPFFits
%
% Run all the PF fits.
%
% The refl option folds the data symmetrically over
% the up-down and down-up configurations.  Sometimes
% this is a cleaner way to look at it, and sometimes not.

%% Clear and close
clear; close all;

%% Some parameters to set
normData = false;
reflectData = false;
corrGuess = false;  

%% 2020-01-31
%
% 7x9 pixel rectangular stimulus, 11043 and 11046
close all;
fitPfsToIncrDecrData('subj','11043','dataDate','20200131','subProject','IncrDecr1','condition','Separation_1', ...
    'norm',normData,'corrGuess',corrGuess,'refl',reflectData);
fitPfsToIncrDecrData('subj','11046','dataDate','20200131','subProject','IncrDecr1','condition','Separation_1', ...
    'norm',normData,'corrGuess',corrGuess,'refl',reflectData);

%% 2021-09-14
%
% 7x9 pixel rectangular stimulus, 11046
% This is a replication of the 2020-01-31 data.
close all;
fitPfsToIncrDecrData('subj','11046','dataDate','20210914','subProject','IncrDecr2','condition','Size_1', ...
    'norm',normData,'corrGuess',corrGuess,'refl',reflectData);

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
% close all;
% fitPfsToIncrDecrData('subj','11046','dataDate','20210914','subProject','IncrDecr2','condition','Size_2', ...
%     'norm',normData,'corrGuess',corrGuess,'refl',reflectData);
% fitPfsToIncrDecrData('subj','11046','dataDate','20210914','subProject','IncrDecr2','condition','Size_2', ...
%     'norm',normData,'corrGuess',corrGuess,'refl',reflectData);

%% 2021-10-18
%
% 5x7 pixel rectangular stimulus, 11046
close all;
fitPfsToIncrDecrData('subj','11046','dataDate','20211018','subProject','IncrDecr3','condition','Size_2', ...
    'norm',normData,'corrGuess',corrGuess,'refl',reflectData);

% 4x8 pixel rectangular stimulus, 11046
% 
% Just one run.  Don't see the incr-decr, and hardly see anything else.
% Too small given gamut limitations.
%
% For some reason I haven't tracked down, the 'refl',false case crashes,
% complaing about violation of the same number of levels in each direction.
% I have not tracked that down.
%
% close all;
% fitPfsToIncrDecrData('subj','11046','dataDate','20211018','subProject','IncrDecr2','condition','Size_3', ...
%     'norm',normData,'corrGuess',corrGuess,'refl',reflectData);
% fitPfsToIncrDecrData('subj','11046','dataDate','20211018','subProject','IncrDecr2','condition','Size_3', ...
%     'norm',normData,'corrGuess',corrGuess,'refl',reflectData);

%% 2021-10-26
%
% 5x7 pixel rectangular stimulus, 11046. 4th quadrant conditions
close all;
fitPfsToIncrDecrData('subj','11046','dataDate','20211026','subProject','IncrDecr4','condition','Size_2', ...
    'norm',normData,'corrGuess',corrGuess,'refl',reflectData);

% 6x8 pixel rectangular stimulus, 11046. 4th quadrant conditions
close all;
fitPfsToIncrDecrData('subj','11046','dataDate','20211026','subProject','IncrDecr4','condition','Size_4', ...
    'norm',normData,'corrGuess',corrGuess,'refl',reflectData);

% 7x9 pixel rectangular stimulus, 11046. 4th quadrant conditions
close all;
fitPfsToIncrDecrData('subj','11046','dataDate','20211026','subProject','IncrDecr4','condition','Size_1', ...
    'norm',normData,'corrGuess',corrGuess,'refl',reflectData);

%% 2021-11-23
%
% 7x9 pixel rectangular stimulus, 11046. 315 degrees, 0, 2, 4, 6, 8 pixel
% separation.  No session normalizing data (0 and 270 deg) collected, which is a shame.
close all;
fitPfsToIncrDecrData('subj','11046','dataDate','20211123','subProject','IncrDecr5','condition','Size1_Sep0', ...
    'norm',normData,'corrGuess',corrGuess,'refl',reflectData);
fitPfsToIncrDecrData('subj','11046','dataDate','20211123','subProject','IncrDecr5','condition','Size1_Sep2', ...
    'norm',normData,'corrGuess',corrGuess,'refl',reflectData);
fitPfsToIncrDecrData('subj','11046','dataDate','20211123','subProject','IncrDecr5','condition','Size1_Sep4', ...
    'norm',normData,'corrGuess',corrGuess,'refl',reflectData);
fitPfsToIncrDecrData('subj','11046','dataDate','20211123','subProject','IncrDecr5','condition','Size1_Sep6', ...
    'norm',normData,'corrGuess',corrGuess,'refl',reflectData);
fitPfsToIncrDecrData('subj','11046','dataDate','20211123','subProject','IncrDecr5','condition','Size1_Sep8', ...
    'norm',normData,'corrGuess',corrGuess,'refl',reflectData);

% 6x8 pixel rectangular stimulus, 11046. 315 degrees, 0 and 4 pixel
% separations.
close all;
fitPfsToIncrDecrData('subj','11046','dataDate','20211123','subProject','IncrDecr5','condition','Size2_Sep0', ...
    'norm',normData,'corrGuess',corrGuess,'refl',reflectData);
fitPfsToIncrDecrData('subj','11046','dataDate','20211123','subProject','IncrDecr5','condition','Size2_Sep4', ...
    'norm',normData,'corrGuess',corrGuess,'refl',reflectData);

%% %% 2021-12-17
%
% 7x9, 11046, separation series. Two runs per direction/separation.
close all;
theConds = {...
    'Dir0_Sep0', ...
    'Dir45_Sep0', ...
    'Dir45_Sep4', ...
    'Dir45_Sep8', ...
    'Dir45_Sep12', ...
    'Dir45_Sep16', ...
    'Dir225_Sep0', ...
    'Dir225_Sep4', ...
    'Dir225_Sep8', ...
    'Dir225_Sep12', ...
    'Dir225_Sep16', ...
    'Dir270_Sep0', ...
    'Dir315_Sep0', ...
    'Dir315_Sep4', ...
    'Dir315_Sep8', ...
    'Dir315_Sep12', ...
    'Dir315_Sep16', ...
    };
for cc = 1:length(theConds)
    fitPfsToIncrDecrData('subj','11046','dataDate','20211217','subProject','IncrDecr6','condition',theConds{cc}, ...
        'norm',normData,'corrGuess',corrGuess,'refl',reflectData);
end

%% 2021-12-17
%
% 7x9, 11002, separation 0, around the ellipse.  Two runs per direction.
close all;
theConds = {...
    'Dir0_Sep0', ...
    'Dir45_Sep0', ...
    'Dir225_Sep0', ...
    'Dir270_Sep0', ...
    'Dir293_Sep0', ...
    'Dir315_Sep0', ...
    'Dir338_Sep0', ...
    };
for cc = 1:length(theConds)
    fitPfsToIncrDecrData('subj','11002','dataDate','20211220','subProject','IncrDecr7_7x9','condition',theConds{cc}, ...
        'norm',normData,'corrGuess',corrGuess,'refl',reflectData);
end

% 6x8, 11002, separation 0. Two runs per direction. Three directions. 
close all;
theConds = {...
    'Dir0_Sep0', ...
    'Dir270_Sep0', ...
    'Dir315_Sep0', ...
    };
for cc = 1:length(theConds)
    fitPfsToIncrDecrData('subj','11002','dataDate','20211220','subProject','IncrDecr7_6x8','condition',theConds{cc}, ...
        'norm',normData,'corrGuess',corrGuess,'refl',reflectData);
end
