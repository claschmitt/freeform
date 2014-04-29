%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% findSpan
%
% find position in the internal knot vector
% INPUT
%   parameter = uk / vl;
%   internalKnotVector = nurbs.knotsU / nurbs.knotsV;
%  optional parameters
%
% OUTPUT
%   span = position of parameter (span) at the knot vector
%
%
% cs, 25.04.2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function span = findSpan(parameter, internalKnotVector)

if parameter == internalKnotVector(1)
    tmpIndex = find(parameter >= internalKnotVector);
    span = tmpIndex(end);
elseif parameter > internalKnotVector(1) && parameter < internalKnotVector(end)
    tmpIndex = find(parameter >= internalKnotVector);
    span = tmpIndex(end);
elseif parameter == internalKnotVector(end)
    tmpIndex = find(parameter == internalKnotVector);
    span = tmpIndex(1) -1;
elseif parameter < internalKnotVector(1)
    disp('Error - findSpan - could not find Span');
    tmpIndex = find(internalKnotVector(1) >= internalKnotVector);
    span = tmpIndex(end);
elseif parameter > internalKnotVector(end)
    disp('Error - findSpan - could not find Span');
    tmpIndex = find(internalKnotVector(end) == internalKnotVector);
    span = tmpIndex(1) -1;
end
