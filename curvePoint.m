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

function C = curvePoint(n,p,U,P,u)

if u >= U(p +1 +n-p +1)
    span = length(U)-p-1 -1;
else
    tmp_span = find(U <= u);
    span = max(tmp_span)-1;
end

N = basisFunction(span,u,p,U);
C = [0.0 0.0 0.0];

C(1,1) = sum (N(0 +1 : p +1) * P(span-p +1 : span-p+p +1,1));
C(1,2) = sum (N(0 +1 : p +1) * P(span-p +1 : span-p+p +1,2));
C(1,3) = sum (N(0 +1 : p +1) * P(span-p +1 : span-p+p +1,3));
%     plot(C(:,1), C(:,2), 'color', 'red');

