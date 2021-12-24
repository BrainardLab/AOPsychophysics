function outAngles = CananicalAngles(inAngles)
% Wrap angles (in degrees) into 0-360.
%
% Syntax:
%    outAngles = CananicalAngles(inAngles)
%
% Description:
%    Wrap the angles so that they 0 <= x < 360.
%
% Inputs:
%    inAngles -     Vector.  Input angles
%
% Outputs:
%    outAngles -    Wrapped version of input.
%
% Optional key/value pairs:
%    None

% History:
%    12/24/21  dhb  Broke out as function.

while (any(inAngles < 0))
    inAngles(inAngles < 0) = inAngles(inAngles < 0) + 360;
end
while (any(inAngles >= 360))
    inAngles(inAngles >= 360) = inAngles(inAngles >= 360) - 360;
end
outAngles = inAngles;