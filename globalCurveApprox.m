%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Global Curve Approximation
%
% compute least squares fitting of a non uniform B-Spline curve with fixed control points
% INPUT
% nurbs = nurbs structure points
% pointsKartesian = observations, measured points[X Y Z]
% Q_ll = kovariance matrix of the observations;
% OUTPUT
% nurbs = nurbs structure
% quality = quality structure
%
% cs, 16.05.2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [nurbs, quality] = globalCurveApprox (nurbs, pointsCartesian, Q_ll)

% compute homologous parameter for measured points
pointsParameter = createPointsParametersCurve(pointsCartesian, 'centripetal');

% compute internal knot vector
nurbs = internalKnots(pointsParameter, nurbs, 'piegl_tiller');

A_full = fillAGM(pointsCartesian, pointsParameter, nurbs);

% Q_ll = sparse(Q_ll);

% reshape pointsCartesian
[numPointsCartesian, dimPoints] = size(pointsCartesian);
pointsRuTMP = reshape(pointsCartesian',numPointsCartesian * dimPoints,1);

% least squares
invN = (A_full' / Q_ll * A_full);
P = invN \ (A_full' / Q_ll * pointsRuTMP);

nurbs.coefs = reshape(P,dimPoints ,nurbs.numberU)';
residuals = A_full * P - pointsRuTMP;

%Quality
quality.Qxx = invN;
%     quality.Qvv = -A * quality.Qxx * A';

quality.sigma0_apost = sqrt((residuals' * residuals) / ...
    (length(pointsRuTMP)-nurbs.numberU)...
    );

quality.residuals = reshape(residuals',dimPoints, numPointsCartesian)';
clear residuals;

quality.sigma_xDach = quality.sigma0_apost .* sqrt(diag(quality.Qxx));


% plot(pointsCartesian(:,1), pointsCartesian(:,2),'o', 'color', 'blue');
% hold on;
% plot(nurbs.coefsU(:,2),nurbs.coefsU(:,3),'-o','color', 'red')
% hold off;


% plot(Q(1:r+1,1), Q(1:r+1,3),'o', 'color', 'blue');
% hold on;
% plot(P_full(:,1),P_full(:,3),'-o','color', 'red')
% hold off;
