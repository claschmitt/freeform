%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basis function
%
% compute parameters for global surface interpolation
% INPUT
% n = control point index function with uk starting with 1 - columns
% m = control point index function with vl starting with 1 - rows
% Q = control Points as a 2 dimensional array with 3 dim coords
%     where each column consist of an 1 dim [x y z]
% OUTPUT
% uk = value in the horizontal knot vector
% vl = value in the vertical knot vector
% 
% cs, 22.05.2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [uk, vl, i_ukvl] = SurfMeshParamsSphere (r,Q,boundingbox)

numSpherePoints = 1000;
if length(Q) > numSpherePoints
    randIndexQ = randi(length(Q),numSpherePoints,1);
else
    randIndexQ = true(length(Q),1);
end

kugel = Kugel_3D(Q(randIndexQ,:));
addpath('Z:\Matlab\ellipsoid_fit');

%% Kugel
G = repmat(kugel.center,r +1,1);
redQ = Q - G;

bbSphere.koord = boundingbox;
bbSphere.koordCentered = bbSphere.koord - G(1:length(bbSphere.koord),:);
for i=1 : length(bbSphere.koord)
    bbSphere.polar(i,1) = atan2(sqrt(bbSphere.koordCentered(i,1)^2 + bbSphere.koordCentered(i,2)^2), bbSphere.koordCentered(i,3)); %vl
    bbSphere.polar(i,2) = atan2(bbSphere.koordCentered(i,2), bbSphere.koordCentered(i,1)); %uk
end

lamda = atan2(sqrt((redQ(:,1) .* redQ(:,1)) + (redQ(:,2) .* redQ(:,2))),redQ(:,3)); %uk
phi   = atan2(redQ(:,2),redQ(:,1)); %vl


%% Norm Lamda & phi
normLamda   = (lamda-min(lamda)) / max(lamda-min(lamda));
uk = normLamda-min(normLamda);

normPhi   = phi / max(phi-min(phi));
vl = normPhi-min(normPhi);
%% delete unique uk vl combinations
[selekted_ukvl i_ukvl i_selekted_ukvl] = unique([uk vl], 'rows');

uk = selekted_ukvl(:,1);
vl = selekted_ukvl(:,2);

