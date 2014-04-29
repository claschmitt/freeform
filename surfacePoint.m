%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute a point on a Nurbs Surface
%
% INPUT
%   nurbs = nurbs structure
%   u = parameter of point in udirection
%   v = parameter of point in v direction
%  optional parameters
%   quality = quality struct of observation
% OUTPUT
%   S = Point on surface
%  optional
%   Qss = Kovariance matrix of the surface point
% 
%
% cs, 25.05.2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [S, varargout] = surfacePoint(nurbs,pointsParameter, varargin)

if isempty(varargin)
    varargout = [];
    controlPointsMatrix = reshapeControlPoints(pointsParameter, nurbs);
    flagQuality = false;
else
    quality = varargin{1};
    [controlPointsMatrix, Qcp] = reshapeControlPoints(pointsParameter, nurbs, quality);
    flagQuality = true;
    
end

AMatrix =  fillAGM([1 1 1],pointsParameter, nurbs);

S = AMatrix * controlPointsMatrix;

if flagQuality
    Qss = A * Qcp * A';
    varargout{1} = Qss;
else
    varargout = [];
end





