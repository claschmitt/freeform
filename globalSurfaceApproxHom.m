%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Global Surface Approximation
%
% compute least squares fitting of a surface with fixed control points
% input
% r = max point index of measured points
% points = mesured points
% p = degree of function 1 with knots U
% q = degree of function 2 with knots V
% n = max index of number of control Points (U-direction)starting with 0
% m = max index of number of control Points (V-direction)starting with 0
% U = knot vektor of function 1 with degree p
% V = knot vektor of function 2 with degree q
%flagImproveParameter = 1 for improvement, 0 no (type: boolean)
%output
% U = knot vektor of function 1 with degree p
% V = knot vektor of function 2 with degree q
% controlP = Control Points of Controle net;
%
% cs, 08.02.2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [nurbs, quality] = globalSurfaceApproxHom (r,pointsCartesian,indexK,nurbs, flagImproveParameter,dataSnooping, boundingbox)

% %% testing 1
% clear all;

%% struct declaration
% quality parameters
quality = struct('sigma0_apost','','sigma_xDach','','residuals','','Qxx','','indexUsedPoints','');
%Nurbs parameters
% nurbs = struct('form','','dimV','','numberU','','dimU','','numberV','','coefs','','orderU','','orderV','','knotsU','','knotsV','');
% DataSnooping
maxIterDataSnooping = 3;
confidenceLevel     = 0.95;
sigma_apri          = 1;
sigma_apri_Koord    = 0.004;
%  Improving observation parameters
maxIterImproveParameters = 1;

%% generate Parameters for measured pointsCartesian

[pointsParameter(:,1), pointsParameter(:,2), quality.indexUsedPoints] = SurfMeshParamsSphere (r,pointsCartesian,boundingbox);
r = length(quality.indexUsedPoints) -1;
pointsCartesian = pointsCartesian(quality.indexUsedPoints,:);

%% DataSnooping

for flagDataSnooping = 1 : maxIterDataSnooping
    nurbs = internalKnots(pointsParameter, nurbs, 'piegl_tiller');
%     %% compute knot vector U
%     r = length(pointsCartesian) -1;
%     if nurbs.knotsU<=0
%         du = (r + 1) / (nurbs.numberU - nurbs.orderU +1);
%         uk_sorted = sort(pointsParameter(:,1));
%         
%         nurbs.knotsU = zeros(nurbs.orderU + nurbs.numberU,1);
%         nurbs.knotsU(nurbs.orderU + nurbs.numberU - nurbs.orderU +1 : nurbs.orderU + nurbs.numberU) = 1;
%         
%         for k=1 : nurbs.numberU-1-nurbs.dimU
%             i = floor(k * du);
%             alpha = k * du - i;
%             %        nurbs.knotsU(nurbs.dimU+k +1) = (1 - alpha) * pointsParameter(:,1)(i-1 +1) + alpha * pointsParameter(:,1)(i +1);
%             nurbs.knotsU(nurbs.dimU+k +1) = (1 - alpha) * uk_sorted(i-1 +1) + alpha * uk_sorted(i +1);
%         end
%     end
%     clear uk_sorted;
%     
%     %% compute knot vector %% compute knot vector U
%     if nurbs.knotsV<=0
%         dv = (r + 1) / (nurbs.numberV - nurbs.orderV +1);
%         vl_sorted = sort(pointsParameter(:,2));
%         
%         nurbs.knotsV = zeros(nurbs.orderV + nurbs.numberV,1);
%         nurbs.knotsV(nurbs.orderV + nurbs.numberV - nurbs.orderV +1 : nurbs.orderV + nurbs.numberV)= 1;
%         
%         for l=1 : nurbs.numberV-1-nurbs.dimV
%             i = floor(l * dv);
%             alpha = l * dv - i;
%             %        nurbs.knotsV(nurbs.dimV+l +1) = (1 - alpha) * pointsParameter(i-1 +1,2) + alpha * pointsParameter(i +1,2);
%             nurbs.knotsV(nurbs.dimV+l +1) = (1 - alpha) * vl_sorted(i-1 +1) + alpha * vl_sorted(i +1);
%         end
%     end
%     
%     clear vl_sorted;
    %% multiplication factor Koordinates range vs. Knot range
    spanCoord = max(pointsCartesian) - min(pointsCartesian);
    spanX     = 1/spanCoord(1);
    spanZ     = 1/spanCoord(3);
    
    %% fill coefficient matrix
    A = fillAGM (pointsCartesian, pointsParameter, nurbs);
    %     tic_fill_A = tic;
    %     %loop to improve parameters
    %     for flagReCalc = 1 : maxIterImproveParameters
    %
    %
    %         A_full = zeros((r +1), (nurbs.numberU-1 +1) * (nurbs.numberV-1 +1));
    %         indexNu_start        = int64(0 +1);
    %         indexNu_end          = int64(nurbs.dimU +1);
    %
    %         for iterAll=0 : r
    %
    %             if flagImproveParameter == true && flagReCalc > 1
    %                 flagBreak_pointsParameter(:,1) = 0;
    %                 flagBreak_vl = 0;
    %
    %                 for iterNewPoint = 1 : 5
    %
    %                     [interpolationPoint, ~] = surfacePoint(nurbs.numberU-1,nurbs.dimU,nurbs.knotsU,nurbs.numberV-1,nurbs.dimV,nurbs.knotsV,nurbs.coefs,pointsParameter(iterAll +1,1),pointsParameter(iterAll +1,2));
    %                     diffPoint = (interpolationPoint - pointsCartesian(iterAll +1,:));
    %
    %                     if abs(diffPoint(1)) > 0.001
    %                         pointsParameter(iterAll +1,1) = pointsParameter(iterAll +1,1) - diffPoint(1)/spanX;
    %                     else
    %                         flagBreak_uk = flagBreak_uk +1;
    %                     end
    %
    %                     if  abs(diffPoint(3)) > 0.001
    %                         pointsParameter(iterAll +1,2) = pointsParameter(iterAll +1,2) + diffPoint(3)/spanZ;
    %                     else
    %                         flagBreak_vl = flagBreak_vl +1;
    %                     end
    %
    %                     if flagBreak_uk == 1 && flagBreak_vl == 1
    %                         break;
    %                     end
    %
    %                 end
    %
    %             end
    %
    %             if pointsParameter(iterAll +1,1) >= nurbs.knotsU(nurbs.numberU-1+1 +1)
    %                 spanU = length(nurbs.knotsU)-nurbs.dimU-1 -1;
    %             else
    %                 tmp_span = find(nurbs.knotsU<=pointsParameter(iterAll +1,1));
    %                 spanU = max(tmp_span)-1;
    %
    %                 clear tmp_span;
    %             end
    %
    %             if pointsParameter(iterAll +1,2) >= nurbs.knotsV(nurbs.numberV-1+1 +1)
    %                 spanV = length(nurbs.knotsV)-nurbs.dimV-1 -1;
    %             else
    %                 tmp_span = find(nurbs.knotsV<=pointsParameter(iterAll +1));
    %                 spanV = max(tmp_span)-1;
    %
    %                 clear tmp_span;
    %             end
    %
    %             Nu = basisFunction(spanU,pointsParameter(iterAll +1,1),nurbs.dimU,nurbs.knotsU);
    %             Nv = basisFunction(spanV,pointsParameter(iterAll +1,2),nurbs.dimV,nurbs.knotsV);
    %
    %             if nurbs.numberU-1<nurbs.dimU
    %                 msgbox('number of controle pointsCartesian must be higher','error');
    %                 break;
    %             end
    %
    %             indexA_row = iterAll +1;
    %
    %             for iterNv = 0 : nurbs.dimV
    %
    %                 indexA_column_start = int64((spanU-nurbs.dimU +1)      +  (nurbs.numberV-1 +1) * (spanV-nurbs.dimV) + iterNv * (nurbs.numberV-1 +1));
    %                 indexA_column_end   = int64((spanU-nurbs.dimU +1) + nurbs.dimU  +  (nurbs.numberV-1 +1) * (spanV-nurbs.dimV) + iterNv * (nurbs.numberV-1 +1));
    %
    %                 A_full(indexA_row,indexA_column_start : indexA_column_end) = Nu(indexNu_start : indexNu_end) * Nv(iterNv +1);
    %
    %
    %             end
    %
    %         end
    %
    %         ElapsedTimeFill_A = toc(tic_fill_A)
    
    %% solve equations
    tic_solve = tic;
    
    %         A = A_full;
    %         clear A_full;
    
    % filling n Matrix
    %         nKlein = zeros((nurbs.numberU-1 +1) * (nurbs.numberV-1 +1),3);
    
    [numPointsCartesian, dimPoints] = size(pointsCartesian);
    l = reshape(pointsCartesian',numPointsCartesian * dimPoints,1);
%     l = zeros(length(pointsCartesian) * 3,1);
%     l(1:3:end) = pointsCartesian(:,1);
%     l(2:3:end) = pointsCartesian(:,2);
%     l(3:3:end) = pointsCartesian(:,3);
    
    nKlein = A' * l;
    
    
    A = sparse(A);
    nKlein = sparse(nKlein);
    
    N = A' * A;
    
    invN = inv(A' * A);
    P = invN * nKlein;
    
    nurbs.coefs      = reshape(P,dimPoints ,nurbs.numberU * nurbs.numberV)';
%     nurbs.coefs(:,1) = P(1:3:end,1);
%     nurbs.coefs(:,2) = P(2:3:end,1);
%     nurbs.coefs(:,3) = P(3:3:end,1);    
    
    ElapsedTimeSolveEquation = toc(tic_solve)
    
    if flagImproveParameter == false
        break
    end
    
    
end
%% return

% compute quality parameters
%     quality = struct('sigma0_apost','','sigma_xDach','','residuals','','Qxx','','indexUsedPoints','');

%residuals
residuals = A * P - l;
clear P;

quality.Qxx = invN;
%     quality.Qvv = -A * quality.Qxx * A';

quality.sigma0_apost = sqrt((residuals' * residuals) / ...
    ((r +1) * 3 - (nurbs.numberV-1 +1) * (nurbs.numberU-1 +1) * 3 )...
    );

quality.residuals = reshape(residuals',dimPoints, numPointsCartesian)';
% quality.residuals = zeros(length(points),3);
% quality.residuals(:,1) = residuals(1:3:end,1);
% quality.residuals(:,2) = residuals(2:3:end,1);
% quality.residuals(:,3) = residuals(3:3:end,1);
clear residuals;

quality.sigma_xDach = quality.sigma0_apost .* sqrt(diag(quality.Qxx));

%Least squares testing
[Ar Ac] = size(A);
degreeFreedom = Ar-Ac;

%chi-square test
%     chiSquare = vartest(quality.sigma0_apost^2,sigma0_apri^2)
chi_quantil    = chi2inv(confidenceLevel,degreeFreedom)
chi_testValue  = degreeFreedom * sigma_apri^2/quality.sigma0_apost^2
chiTest_result = chi_testValue < chi_quantil;

%t-test
t_quantil = tinv(confidenceLevel,degreeFreedom);
t_testValue  = abs(quality.residuals/sigma_apri_Koord);
errorIndex = t_testValue < t_quantil;
errorIndex = logical(errorIndex(:,1) .* errorIndex(:,2) .* errorIndex(:,3));



if dataSnooping  && flagDataSnooping < maxIterDataSnooping && (sum(errorIndex) < length(errorIndex)-10)
    disp(['Iteration data snooping: ' num2str(flagDataSnooping)]);
    
    points = points(errorIndex,:);
    r = length(points) -1;
    pointsParameter(:,1) = pointsParameter(errorIndex,1);
    pointsParameter(:,2) = pointsParameter(errorIndex,2);
    nurbs.knotsU = 0;
    nurbs.knotsV = 0;
    quality.indexUsedPoints = quality.indexUsedPoints (errorIndex,1);
else
    %         break;
end
%% plot results

% plot3(points(:,1), points(:,2), points(:,3),'o', 'color', 'blue');

% hold on;
% plot3(nurbs.coefs(:,1),nurbs.coefs(:,2),nurbs.coefs(:,3),'-o','color', 'red')
% hold off;
