
%test globalSurfApprox

%% load Kappe
clear all;
close all;
clc;
%Load Data
    scan = ImportXYZI('Z:\12_Forschung\02_Projekte\02_Internal\FreeForm\daten\2_K2_u_cut2.pts', 'headerline', 1);

% G2 - K2
     %ALL
    scan.boundingBox = [105.010 106.530 110.6169 ;...
                        110.060 106.530 110.6169 ;...
                        107.480 100.000 114.500 ;...
                        105.010 106.530 110.6169];
    
%Konvexe Hülle

    insideXY = inpolygon(scan.data{1}(:,1), scan.data{1}(:,2), scan.boundingBox(:,1),  scan.boundingBox(:,2)); 
    insideZ = (scan.data{1}(:,3) > min( scan.boundingBox(:,3))) & (scan.data{1}(:,3) < max( scan.boundingBox(:,3)));
    inside = logical(insideXY .* insideZ);
%     clear insideXY;
%     clear insideZ;
        
    pointsInside = scan.data{1}(inside,1:3);
%     clear inside;

   [indexK , volume] = convhull(pointsInside(:,1), pointsInside(:,2)); %2D

% points = load(cloudpath);
    roundFactor = 350;
        tmpPoints = floor(pointsInside(:,1:3) * roundFactor)/roundFactor;
%     clear scan.data;
    tmpPoints = unique(tmpPoints,'rows','stable');
    points = vertcat(tmpPoints, pointsInside(indexK,1:3));
%     clear tmpPoints pointsInside;
             
    meanPoint = mean(points);
    meanPoints = repmat(meanPoint,length(points),1);
    points = points - meanPoints;
%     clear meanPoints;
    

% for inm=15 : 20
    inm = 10
%Definitions
    nurbs = struct('form'   ,'B-Spline-Surface',...
                   'dimV'   ,         3,... %degree of function in U direction
                   'dimU'   ,         1,... %degree of function in V direction
                   'numberU',        11,... % number of controle points in u direction
                   'numberV',        11,... % number of controle points in v direction
                   'coefs'  ,        '',...
                   'orderU' ,         2,... %order of function in U direction
                   'orderV' ,         4,... %order of function in V direction
                   'knotsU' ,         0,...
                   'knotsV' ,         0);
    % Degree 
    p = 2; %degree of function in U direction
    q = 2; %degree of function in V direction

    % num measured points
    r = length(points) -1; %index of the measured points in u direction
    s = r; %index of the measured points in v direction

    % num control points 
    n = inm; %max index of controle points in u direction
    m = inm; %max index of controle points in v direction
   
%% Global Surface Approximation Homogeneous
% -----
% testing approx

% -----
numPoints4calculation = 1000000;
% numPoints4calculation = length(points);

if length(points) > numPoints4calculation
    randIndexPoints = randi((r +1), numPoints4calculation,1);
    points = points(randIndexPoints,:);
    r = length(points) -1;
end

[nurbs, quality] = globalSurfaceApproxHom (r,points,indexK,nurbs,false,false, scan.boundingBox);

%nurbs control Points + barycentric center
nurbs.coefs = nurbs.coefs + repmat(meanPoint,length(nurbs.coefs),1);

% -----
% testing plot
% -----
figureSurface = figure('name', 'Surface plot');
hold on;
resolution = 50;
[intPoint, intPoint_normal] = plotSurfaceHom(nurbs,resolution, figureSurface, [], points(quality.indexUsedPoints,:), quality);

MatLabNurbs = WrapperNu2MaNu(nurbs);

Name = [num2str(p) '_' num2str(n +1)];
%Nurbs export
igesout(MatLabNurbs,[Name '_cut_NURBS-NURBS' ]);
% end
