function [pCorrected] = CorrectForGuessing(pHit,pFA)
% Correct for guessing
%
% Synopsis:
%
% Description:
%    Apply standard correction for guessing formula.
%
%    False alarm rate is corrected to 0, and anything less
%    than that also is forced to 0, since negative doesn't
%    make sense.
%
% Inputs:
%    pHit -         Vector of hit rates.
%    pFA -          False alarm rate to use in the correction
%
% Outputs:
%   pCorrected -    Vector. pHit corrected for guessing.

% History:
%   12/23/21 dhb    Pulled out into its very own function. 
%                   Added avoidance of negative return values.

% Avoid the weird world of negative proportion correct.
pHit(pHit < pFA) = pFA;

% Standard correct for guessing formula.
pCorrected = (pHit-pFA)/(1-pFA);

end