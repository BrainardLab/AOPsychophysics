% RunAllPFFits
%
% Run the full set of initial and aggregated PF Fits.

% Initial fits.  This is set so that it doesn't correct
% for guessing, and the main point of running it is to massage the data to
% allow the aggregated fitting calls below to operate.  Correction for
% guessing and slope constraining is done there, across multiple runs with
% different stimulus parameters, as controlled by those routines.
clear; close all; RunAllInitialPFFits;

% Aggregated fits.  Correction for guessing is done within session
% separatly for each stimulus size, and psychometric function slope is
% constrained across stimuli of the same size.  This can be fairly finely
% controlled if desired, but this approach seems as sensible as anything.
%
% Note the PF fits can get stuck, so that grid search is important.
% Sometimes the parameters for the grid search need to be tweaked by hand
% after an initial run, so that fits converge sensibly.  Worth eyeballing
% the fit PFs after these are run.
%
% I am not sure how Palamedes weights fits across different numbers of
% trials in different directions, but not inclined to delve too far down
% that rabit hole.
clear; close all; AggregatePFFits_5by7_11046;
clear; close all; AggregatePFFits_6by8_11002;
clear; close all; AggregatePFFits_6by8_11046;
clear; close all; AggregatePFFits_7by9_11002;
clear; close all; AggregatePFFits_7by9_11043;
clear; close all; AggregatePFFits_7by9_11046;