%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Distance3D
%
% compute the euclidian distance of two three dimensional points
% INPUT
% p1 = point 1 as a array [x y z]
% p2 = point 2 as a array [x y z]
% % OUTPUT
% d =
%  
% cs, 22.05.2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [d] = distPoint2Point(p1, p2)

[rp1 cp1] = size(p1);
[rp2 cp2] = size(p2);


if rp1 ==1 && rp2==1 &&  (cp1==2 || cp2==2 )
    d = sqrt( (p2(1)-p1(1))*(p2(1)-p1(1)) + (p2(2)-p1(2))*(p2(2)-p1(2)));
elseif rp1 ==1 && rp2==1 &&  cp1==3 && cp2==3
    d = sqrt( (p2(1)-p1(1))*(p2(1)-p1(1)) + (p2(2)-p1(2))*(p2(2)-p1(2)) + (p2(3)-p1(3))*(p2(3)-p1(3)) );
else
    error('Matrix dimension not valid');
end
