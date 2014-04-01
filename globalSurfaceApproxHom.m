%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Global Surface Approximation
%
% compute least squares fitting of a surface with fixed control points
% input
% r = max point index of measured points 
% Q = mesured points
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


function [nurbs, quality] = globalSurfaceApproxHom (r,points,indexK,q,p,n,m,uk,vl,U,V, flagImproveParameter,dataSnooping, boundingbox)

% %% testing 1
% clear all;

%% struct declaration
% quality parameters
quality = struct('sigma0_apost','','sigma_xDach','','residuals','','Qxx','','indexUsedPoints','');
%Nurbs parameters
nurbs = struct('form','','dimV','','numberU','','dimU','','numberV','','coefs','','orderU','','orderV','','knotsU','','knotsV','');
% DataSnooping
maxIterDataSnooping = 3;
confidenceLevel     = 0.95;
sigma_apri          = 1;
sigma_apri_Koord    = 0.004;
%  Improving observation parameters   
maxIterImproveParameters = 1;

%% generate Parameters for measured points

[uk, vl, quality.indexUsedPoints] = SurfMeshParamsSphere (r,points,boundingbox);
r = length(quality.indexUsedPoints) -1;
Q = points(quality.indexUsedPoints,:);
clear points;

%% DataSnooping

for flagDataSnooping = 1 : maxIterDataSnooping
    %% compute knot vector U
    r = length(Q) -1;
    if U<=0
        du = (r + 1) / (n - p + 1);
        uk_sorted = sort(uk);
        
        U = zeros((p+1)*2+n-p,1);
        U (p+1+n-p +1:(p+1)*2+n-p)= 1;
        
        for k=1 : n-p
            i = floor(k * du);
            alpha = k * du - i;
            %        U(p+k +1) = (1 - alpha) * uk(i-1 +1) + alpha * uk(i +1);
            U(p+k +1) = (1 - alpha) * uk_sorted(i-1 +1) + alpha * uk_sorted(i +1);
        end
    end
    clear uk_sorted;
    
    %% compute knot vector %% compute knot vector U
    if V<=0
        dv = (r + 1) / (m - q + 1);
        vl_sorted = sort(vl);
        
        V= zeros((q+1)*2+m-q,1);
        V (q+1+m-q +1:(q+1)*2+m-q)= 1;
        
        for l=1 : m-q
            i = floor(l * dv);
            alpha = l * dv - i;
            %        V(q+l +1) = (1 - alpha) * vl(i-1 +1) + alpha * vl(i +1);
            V(q+l +1) = (1 - alpha) * vl_sorted(i-1 +1) + alpha * vl_sorted(i +1);
        end
    end
    
    clear vl_sorted;
    %% multiplication factor Koordinates range vs. Knot range
    spanCoord = max(Q) - min(Q);
    spanX     = 1/spanCoord(1);
    spanZ     = 1/spanCoord(3);
    
    %loop to improve parameters
    for flagReCalc = 1 : maxIterImproveParameters
        %% fill coefficient matrix in U direction without the first and last points
        tic_fill_A = tic;
        
        A_full = zeros((r +1), (n +1) * (m +1));
        indexNu_start        = int64(0 +1);
        indexNu_end          = int64(p +1);
        
        for iterAll=0 : r
            
            if flagImproveParameter == true && flagReCalc > 1
                flagBreak_uk = 0;
                flagBreak_vl = 0;
                
                for iterNewPoint = 1 : 5
                    
                    [interpolationPoint, ~] = surfacePoint(n,p,U,m,q,V,P,uk(iterAll +1),vl(iterAll +1));
                    diffPoint = (interpolationPoint - Q(iterAll +1,:));
                    
                    if abs(diffPoint(1)) > 0.001
                        uk(iterAll +1) = uk(iterAll +1) - diffPoint(1)/spanX;
                    else
                        flagBreak_uk = flagBreak_uk +1;
                    end
                    
                    if  abs(diffPoint(3)) > 0.001
                        vl(iterAll +1) = vl(iterAll +1) + diffPoint(3)/spanZ;
                    else
                        flagBreak_vl = flagBreak_vl +1;
                    end
                    
                    if flagBreak_uk == 1 && flagBreak_vl == 1
                        break;
                    end
                    
                end
                
            end
            
            if uk(iterAll +1) >= U(n+1 +1)
                spanU = length(U)-p-1 -1;
            else
                tmp_span = find(U<=uk(iterAll +1));
                spanU = max(tmp_span)-1;
                
                clear tmp_span;
            end
            
            if vl(iterAll +1) >= V(m+1 +1)
                spanV = length(V)-q-1 -1;
            else
                tmp_span = find(V<=vl(iterAll +1));
                spanV = max(tmp_span)-1;
                
                clear tmp_span;
            end
            
            Nu = basisFunction(spanU,uk(iterAll +1),p,U);
            Nv = basisFunction(spanV,vl(iterAll +1),q,V);
            
            if n<p
                msgbox('number of controle points must be higher','error');
                break;
            end
            
            indexA_row = iterAll +1;
            
            for iterNv = 0 : q
                
                indexA_column_start = int64((spanU-p +1)      +  (m +1) * (spanV-q) + iterNv * (m +1));
                indexA_column_end   = int64((spanU-p +1) + p  +  (m +1) * (spanV-q) + iterNv * (m +1));
                
                A_full(indexA_row,indexA_column_start : indexA_column_end) = Nu(indexNu_start : indexNu_end) * Nv(iterNv +1);

                
            end
            
        end
        
        ElapsedTimeFill_A = toc(tic_fill_A)
        
        %% solve equations
        tic_solve = tic;
        
        A = A_full;
        clear A_full;
        
        % filling n Matrix
%         nKlein = zeros((n +1) * (m +1),3);
        nKlein = A' * Q;

        
        A = sparse(A);
        nKlein = sparse(nKlein);

        N = A' * A;

        invN = inv(A' * A);
        P = invN * nKlein;


        ElapsedTimeSolveEquation = toc(tic_solve)
        
        if flagImproveParameter == false
            break
        end
        
        
    end
    %% return
    
    % compute quality parameters
%     quality = struct('sigma0_apost','','sigma_xDach','','residuals','','Qxx','','indexUsedPoints','');
    
    %residuals
    quality.residuals = A * P - Q;
    
    quality.Qxx = invN;
%     quality.Qvv = -A * quality.Qxx * A';
    
    quality.sigma0_apost = sqrt(([quality.residuals(:,1); quality.residuals(:,2); quality.residuals(:,3)]' * ...
        [quality.residuals(:,1); quality.residuals(:,2); quality.residuals(:,3)]) / ...
        ((r +1) * 3 - (m +1) * (n +1) * 3 )...
        );
    
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

        Q = Q(errorIndex,:);
        r = length(Q) -1;
        uk = uk(errorIndex);
        vl = vl(errorIndex);
        U = 0;
        V = 0;
        quality.indexUsedPoints = quality.indexUsedPoints (errorIndex,1);
    else
        break;
    end
    
end
% Nurbs Values
% nurbs = struct('form','','dimV','','numberU','','dimU','','numberV','','coefs','','orderU','','orderV','','knotsU','','knotsV','');
nurbs.form = 'B-Spline-Surface';
nurbs.dimV  = q;
nurbs.numberU = n +1;
nurbs.dimU  = p;
nurbs.numberV = m +1;
nurbs.coefs  = P';
nurbs.orderU  = p +1;
nurbs.orderV  = q +1;
nurbs.knotsU  = U;
nurbs.knotsV  = V;

%% plot results

% P_full = zeros(n+1,3);
% P_full(1,:) = Q(1,:);
% P_full(n+1,:) = Q(r+1,:);
% P_full(2:n,1) = P(:,1); 
% P_full(2:n,2) = P(:,2);
% P_full(2:n,3) = P(:,3);

% plot3(Q(:,1), Q(:,2), Q(:,3),'o', 'color', 'blue');

% hold on;
% plot3(P(:,1),P(:,2),P(:,3),'-o','color', 'red')
% hold off;
