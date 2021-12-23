function outVec = MatchEntriesToTolerance(inVec,tolerance)
% Match up entries in the input that are within tolerance
%
% Syntax:
%    outVec = MatchEntriesToTolerance(inVec,tolerance)
%
% Description:
%    Sometimes we want to match up entries in a vector that are within
%    tolerance of each other, treating them as the same number.  This
%    routine goes through the input and changes entries that are within
%    tolerance to be the same as each other.  Which one gets adjusted is
%    not specified; that is, if you're using this routine you shouldn't
%    care about that.
%
%    An application of this is that in some experiments where we specifiy
%    vector direction of stimuli, a change in angle occurs somewhere
%    because of rounding.  This routine can be used to patch that up.
%
% Inputs:
%    inVec       - The input vector
%    tolerance   - Match up to +/- this tolerance
%
% Outputs:
%    outVec      - The patched up output vector.
%
% Optional key/value pairs
%    None
%
% See also:
%

% History:
%    12/23/21  dhb  Wrote it

% Examples:
%{
    inVec = [1 5 7 8 10];
    tolerance = 1;
    outVec = MatchEntriesToTolerance(inVec,1)
%}

outVec = inVec;
for ii = 2:length(inVec)
    for jj = 1:ii-1
    	if (abs(inVec(ii)-inVec(jj)) <= tolerance)
            outVec(ii) = outVec(jj);
        end
    end
end