function outAngles = ReflectAnglesAround45(inAngles)
% Reflect angles to canonical side of 45-deg line
%
% Synopsis:
%    outAngles = ReflectAnglesAround45(inAngles)
% 
% Description:
%    In the two-spot experiments, stimuli that are symmetric around the
%    45-deg line are the same except for an up-down/down-up reversal.
%
%    This routine relfects angles 45 < x < 225 to the other side of the 45
%    degree line.
%
%    Angles in degrees.
%
% Inputs:
%
% Outputs:
%
% Optional key/value pairs
%    None
%
% See also: CanonicalAngles, MatchEntriesToTolerance.

% History:
%   12/25/21  dhb  Wrote it.

outAngles = inAngles;
for aa = 1:length(inAngles)
    if (inAngles(aa) > 45 && inAngles(aa) < 225)
        delta = inAngles(aa) - 45;
        outAngles(aa) = 45 - delta;
    end
end
outAngles = CanonicalAngles(outAngles);