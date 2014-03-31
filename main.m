%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mainFIG
%
%
% input
%
%
% out
%test
%
%
% cs, 19.03.2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;
close all;

% add path for other functions

%load scanfile
load('Z:\12_Forschung\02_Projekte\02_Internal\FreeForm\daten\Nr14183 ultrahigh I637 Profiler.part10.xyz.asc.mat');
xyz = scan.data{1,1}(:,1:3);
xyz(:,2) = xyz(:,2) - abs(min(scan.boundingBox(:,1)));
clear scan;

% % calculate observation values
% % hzVD = XYZ2angleDistance (xyz);
% load('Nr14183 ultrahigh I637 Profiler.part10.xyz.asc - Profil1 hzVD.mat');
% 
% %calculate incidence angle
% % iA = IncidenceAngle([0 0], xyz(:,2:3), 1);
% load('Nr14183 ultrahigh I637 Profiler.part10.xyz.asc - Profil1 IncidenceAngle.mat');

% reduce Data
boundaryMin = -4;
boundaryMax = +4;
indexPointsBoundary = find(xyz(:,2)> boundaryMin & xyz(:,2) < boundaryMax);
xyz  = xyz (indexPointsBoundary,:);

% iA   = iA  (indexPointsBoundary,:);
% hzVD = hzVD(indexPointsBoundary,:);

% numPoints = 'all'; % all Points
% numPoints = 8000;
% if not(strcmp('all', numPoints ))
%     xyz  = xyz (1:numPoints,:);
%     iA   = iA  (1:numPoints,:);
%     hzVD = hzVD(1:numPoints,:);
% end


%% Uncertainty by Alkhatib 2009, Lisabon
% Parameter in [m] or [rad]
% parameter = struct (...
%     'ad',                    0.0000, ...
%     'iA',                        iA, ...
%     'zenitI',                0.0000, ...
%     'zenitR',                0.0000, ...
%     'sigma_d',               0.0010, ... %0.0030, ...
%     'sigma_ad',              0.0000, ...
%     'sigma_ppm',               2E-5, ...
%     'sigma_iA',          iA * 0.005, ... %?
%     'sigma_zenit',   0.005 * pi/200, ...
%     'sigma_zenitI',  0.000 * pi/200, ...
%     'sigma_zenitR',  0.020 * pi/200 );
% 
% %calculate kofactor matrix for xyz values
% Q_YZ = VfTLSProfile (hzVD,parameter);

%% plotting
%plot accuracies
% mainDiag_Q_YZ = sqrt(diag(Q_YZ));
% figure
% title('Accuracy TLS Datapoints');
% hold on;
% 
% plot(xyz(:,2),mainDiag_Q_YZ(1:2:end),'r');
% [AX, H1, H2] = plotyy(xyz(:,2),mainDiag_Q_YZ(2:2:end),xyz(:,2),abs(iA) * 200/pi,'plot');
% % plot(xyz(:,2),abs(iA)/1000,'m');
% 
% hold off;
% legend ('y coordinate', 'z coordinate', 'incidence angle');
% set(get(AX(1), 'xlabel'), 'String', 'Position at the bridge [m]');
% set(get(AX(1), 'ylabel'), 'String', 'Accuracy [mm]');
% set(get(AX(2), 'ylabel'), 'String', 'Incidence angle [gon]');

%% B-Spline fitting

% extract linear trend
trendParams = polyfit(xyz(:,2), xyz(:,3),1);
xyzRed(:,3) = (xyz(:,3) - polyval(trendParams,xyz(:,2))) * 1000;

% make positive
minY = min(xyz(:,2));
xyzRed(:,2) = xyz(:,2) - minY;

%         number of points
r = length(xyzRed);
r = r-1;
%         number of points in orthogonal to r
s = 0;
%         degree of B-Spline
p = 3;
%         max index of internal knots
n = 20;

[uk vl] = SurfMeshParamsUnstrukt (r,s,[xyzRed(:,2) zeros(r +1,1) xyzRed(:,3)]);

%         if U<=0
% compute knot vector U
du = (r + 1) / (n - p + 1);

U = zeros((p+1)*2+n-p,1);
U (p+1+n-p +1:(p+1)*2+n-p)= 1;

for k=1 : n-p
    i = floor(k * du);
    alpha = k * du - i;
    U(p+k +1) = (1 - alpha) * uk(i-1 +1) + alpha * uk(i +1);
end


% [U_VK,P_VK] = globalCurveApprox (r,[xyzRed(:,2) zeros(r +1,1) xyzRed(:,3)],p,n,uk,U,Q_YZ);
% 
% %re translate & re trend
% P_VK(:,1) = P_VK(:,1) + minY;
% % P_VK(:,3) = P_VK(:,3) + polyval(trendParams,P_VK(:,1));
% 
% nurbs_VK.form = 'B-Spline';
% nurbs_VK.dimU  = p;
% nurbs_VK.numberU = n +1;
% nurbs_VK.coefs  = P_VK';
% nurbs_VK.orderU  = p+1;
% nurbs_VK.knotsU  = U_VK;
% 
% 
% % Plotting profile
% figureCurve_VK = figure;
% plot(xyzRed(:,2) + minY, xyzRed(:,3),'.g')
% hold on;
% subPoints_VK = plotCurve(nurbs_VK,1000,figureCurve_VK);
% hold off;
% title('B-Spline with VK information')
% legend('Measured points','B-spline', 'Control net');
% xlabel('x-coordinate [m]');
% ylabel('red. z-coordinate [mm]');
% %   fileNameFigure = ['Y:\Forschung\AusIng\institut_fuer_massivbau\rethen\20130701_Messung 3\Zeitreihen\' fileNameTLS '_Profile'   num2str(counterProfile) '_.fig'];
% %   saveas(figureCurve,fileNameFigure, 'fig');
% %   close(figureCurve);


[U,P] = globalCurveApprox (r,[xyzRed(:,2) zeros(r +1,1) xyzRed(:,3)],p,n,uk,U,eye((r+1)*2));
%re translate & re trend
P(:,1) = P(:,1) + minY;
% P(:,3) = P(:,3) + polyval(trendParams,P(:,1));

nurbs.form = 'B-Spline';
nurbs.dimU  = p;
nurbs.numberU = n +1;
nurbs.coefs  = P';
nurbs.orderU  = p+1;
nurbs.knotsU  = U;


% Plotting profile
figureCurve = figure;
plot(xyzRed(:,2) + minY, xyzRed(:,3),'.g')
hold on;
subPoints = plotCurve(nurbs,1000,figureCurve);
hold off;
title('B-Spline with EYE information')
legend('Measured points','B-spline', 'Control net');
xlabel('x-coordinate [m]');
ylabel('red. z-coordinate [mm]');
xlim([boundaryMin boundaryMax]);
%   fileNameFigure = ['Y:\Forschung\AusIng\institut_fuer_massivbau\rethen\20130701_Messung 3\Zeitreihen\' fileNameTLS '_Profile'   num2str(counterProfile) '_.fig'];
%   saveas(figureCurve,fileNameFigure, 'fig');
%   close(figureCurve);

% % Plot B-Splines VK + EYE and
% figureBuVK_X = figure;
% plot(subPoints_VK(:,1),subPoints_VK(:,3), '-b');
% hold on;
% plot(subPoints(:,1),subPoints(:,3), '-r');
% title('B-Spline with and without VK information')
% legend('B-Spline VK','B-Spline EYE');
% xlabel('x-coordinate [m]');
% ylabel('red. z-coordinate [mm]');
% hold off;
% 
% % differenc of both profiles with VK and Eye
% diff_VK_EYE = subPoints - subPoints_VK;
% figureDiff_X = figure;
% plot((subPoints(:,1) + subPoints_VK(:,1)) ./2 , diff_VK_EYE(:,1), '+b');
% hold on;
% % figureDiff_Z = figure;
% plot((subPoints(:,1) + subPoints_VK(:,1)) ./2 , diff_VK_EYE(:,3), '+r');
% title('Difference B-Spline with and without VK information')
% legend('Diff X-coordinate','Diff Z-coordinate');
% xlabel('x-coordinate [m]');
% ylabel('red. z-coordinate [mm]');
% hold off;







