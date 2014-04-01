%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute a point on a Nurbs Surface
%
% INPUT
% n = max index of control points in u direction
% p = degree of function in u direction
% U = knot vector in u direction
% m = max index of control points in v direction
% q = degree of function in v direction
% V = knotvector in v direction
% P = controlpoints
% u = parameter of point in udirection
% v = parameter of point in v direction
% OUTPUT
% S = Point on surface
% dersS = value of derivate of the surface point
% 
% cs, 25.05.2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [S] = surfacePoint(n,p,U,m,q,V,P,u,v)

if u >= U(p +1 +n-p +1)
    spanu = length(U)-p-1 -1;
else
    tmp_spanu = find(U <= u);
    spanu = max(tmp_spanu)-1;
end
Nu = basisFunction(spanu,u,p,U);


if v >= V(q +1 + m-q +1)
    spanV = length(V)-q-1 -1;
else
    tmp_spanV = find(V <= v);
    spanV = max(tmp_spanV)-1;
end
Nv = basisFunction(spanV,v,q,V);
    
if sum(isnan(Nu))== 0 && sum(isnan(Nv))== 0
    
   
    tmpu = zeros(q +1,3);
    derstmpu = zeros(q +1,3);
    
    for l=0 : q
    % calculation is first the u direction
        
        indexv = (m +1) * (spanV-q) + l * (m +1);
        indexu = (spanu-p);
       
        tmpu(l +1,:) = Nu * P(indexv + indexu +1 : indexv + indexu +p +1, :);        

    end
    
    S = Nv * tmpu;

   
else
    S = [];
end




