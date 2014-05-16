%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute a point on the Spline
%
% INPUT
% n = max index control points
% p = degree
% U = knot vektor
% P = ControlPoints
% u = new Point on the NURBS Value between 0 and 1
% OUTPUT
% C = Point on curve
%
% cs, 25.05.2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function curvePoint = curvePoint(pointsParameter, nurbs)

[r, dimPoints] = size(nurbs.coefs);
spanU = findSpan(pointsParameter,nurbs.knotsU);

coefsReshaped = reshape(nurbs.coefs(spanU-nurbs.orderU +1 : spanU-nurbs.orderU + nurbs.orderU,:)',nurbs.orderU * dimPoints,1);
AMatrix = fillAGM (zeros(1,dimPoints), pointsParameter, nurbs);

curvePoint = (AMatrix * coefsReshaped)'; 

