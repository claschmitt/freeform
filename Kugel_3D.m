%**********************************************************************************************
% Ausgleichung im GH-Modell (streng)
% Ausgleichender Kugel 3D
%
% Observations   : yi, xi, zi
% Unknowns       : Ym, Xm, Zm, R
% 
% INPUT
% points3D = [x y z]
% varargin = {'Sigma0' Wert ;
%             'MaxD' Wert; (Abbruchkriterium Maximales delata der abgezogenen Beobachtungen)
%             'Iterations' Wert
%             'Datasnooping' number of thrown out Points on each iteration
%             'ReportPath'
%             'Q_ll' Kovarianzmatrix of the observations a priori
%             'Plot': true/fals}
% 
% OUTPUT 
% kugel = struct('center',[Xm_dach Ym_dach Zm_dach],...
%                'sigma_center', [m_Xm m_Ym m_Zm]  ,...
%                'radius', R_dach                  ,...
%                'sigma_radius', m_R               ,...
%                'residuals',v_dach                ,...
%                'Qxx',Q_xx                        ,...
%                'sigma0_apost',sigma_0_a_post     ,...
%                'iterations', iteration           ,...
%                'balancedCoord', [x + v_0_x + gravityC(1), y + v_0_y + gravityC(2), z + v_0_z + gravityC(3)]);
%
% Author         : Claudius Schmitt
% Version        : 08.02.2013
% Last changes   : 
%                  
%**********************************************************************************************

function [kugel] = Kugel_3D(point3D, varargin)

%% Test
%clear all; clc; close all;
%addpath('Z:\Matlab\NURBS');
%   scan = ImportXYZI('Y:\Forschung\AusIng\institut_fuer_massivbau\neugotischekirche_hildesheim\Gewoelbemodellierung\Punktwolken\Gu2_west_u.pts','headerline',1);
%randIndexQ = floor(rand(3000,1) * length(scan.data{1}));
%point3D = scan.data{1}(randIndexQ,:);
% load('points_flaeche2Grades_test.mat')
%varargin = {'Plot' 1 ;...
%            'Datasnooping', 20;...
%            'ReportPath' 'Y:\Forschung\AusIng\institut_fuer_massivbau\neugotischekirche_hildesheim\Gewoelbemodellierung\Punktwolken\Kugel3D-Gu2_west_u.out';...
%            'sigma0' 1};



%%
tic_sphere = tic;

%default values overwriten when denoted on function call
sigma_0_a_pri = 0.0005;
EPSILON = 10e-3;
maxIterations = 20;
dataSnooping = 1;
plot = false;
confidenceLevel = 0.975; %Zweiseitig 0.95 einseitig

if isempty(varargin)
    
else
    
    indexSigma0 = strcmp('sigma0', varargin);
    if isempty(find(indexSigma0,1))
        sigma_0_a_pri = 1.0;
    else
        sigma_0_a_pri = varargin{find(indexSigma0,1), 2};
    end

    indexAbbr   = strcmp('MaxD', varargin);
    if isempty(find(indexAbbr,1))
        EPSILON = 10e-3;
    else
        EPSILON = varargin{find(indexAbbr,1), 2};
    end

    indexIterations   = strcmp('Iterations', varargin);
    if isempty(find(indexIterations,1))
        maxIterations = 10;
    else
        maxIterations = varargin{find(indexIterations,1), 2};
    end
    
    indexDatasnooping   = strcmp('Datasnooping', varargin);
    if isempty(find(indexDatasnooping,1))
        dataSnooping = 1;
    else
        dataSnooping = varargin{find(indexDatasnooping,1), 2};
    end
    
    indexReportPath   = strcmp('ReportPath', varargin);
    if isempty(find(indexReportPath,1))
        reportPath = 'Kugel_3D.out';
    else
        reportPath = varargin{find(indexReportPath,1), 2};
    end
    
    indexQ_ll   = strcmp('Q_ll', varargin);
    if isempty(find(indexQ_ll,1))
        Q_ll = [];
    else
        Q_ll = varargin{find(indexQ_ll,1), 2} * sigma_0_a_pri^2;
    end
    
    indexPlot   = strcmp('Plot', varargin);
    if isempty(find(indexPlot,1))
        plot = false;
    else
        plot = varargin{find(indexPlot,1) +1};
    end
    
end

iSnoop = 1;
while iSnoop <= dataSnooping
    iteration = 0;
    max_delta_x = 10e10;
    
    gravityC = mean(point3D(:,1:3));

    x = point3D(:,1) - gravityC(1);
    y = point3D(:,2) - gravityC(2);
    z = point3D(:,3) - gravityC(3);
%     clear point3D;

    l = vertcat(x,y,z);

    % Anzahl der Messwerte ermitteln
    n  = length(l);
    nx = n / 3;
    % Startwerte fuer die Unbekannten
    Ym_0 = 2*sum(y)/n;
    Xm_0 = 2*sum(x)/n;
    Zm_0 = 2*sum(z)/n;
    R_0  = 7.7;
%     R_0  = 9;

    % Startwerte fuer die Verbesserungen

    % Vektor v_0
    v_0   = zeros(n,1);
    v_0_x = zeros(nx,1);
    v_0_y = zeros(nx,1);
    v_0_z = zeros(nx,1);

    % Kofaktorenmatrix der Beobachtungen
    if exist('Q_ll','var')
        if isempty(Q_ll) 
            Q_ll = eye(n,n) * sigma_0_a_pri^2;
        end
    else
        Q_ll = eye(n,n) * sigma_0_a_pri^2;
    end

    % --> Start iteration least squares

    while max_delta_x > EPSILON

        repmat_Xm_0 = repmat(Xm_0,nx,1);
        repmat_Ym_0 = repmat(Ym_0,nx,1);
        repmat_Zm_0 = repmat(Zm_0,nx,1);

        % Matrix B
        dpsi_dvx = zeros(1,nx);
        dpsi_dvy = zeros(1,nx);
        dpsi_dvz = zeros(1,nx);

        dpsi_dvx(1,:) = 2*((x(:,1) + v_0_x) - repmat_Xm_0);
        dpsi_dvy(1,:) = 2*((y(:,1) + v_0_y) - repmat_Ym_0);
        dpsi_dvz(1,:) = 2*((z(:,1) + v_0_z) - repmat_Zm_0);

        B1 = diag(dpsi_dvx);
        B2 = diag(dpsi_dvy);
        B3 = diag(dpsi_dvz);

        B = horzcat(B1,B2,B3);

        % Matrix A
        A = zeros(nx,4);
        A(:,1) = -2*((x(:,1) + v_0_x) - repmat_Xm_0);
        A(:,2) = -2*((y(:,1) + v_0_y) - repmat_Ym_0);
        A(:,3) = -2*((z(:,1) + v_0_z) - repmat_Zm_0);
        A(:,4) = repmat(-2*R_0,nx,1);

        % Bedingungen als Funktion der NW fuer die Unbekannten und Verbesserungen
        b_0 = zeros(nx,1);

        % Bedingungen als Funktion der NW fuer die Unbekannten und Verbesserungen
        b_0(:,1) =  (((x(:,1) + v_0_x) - repmat_Xm_0).^2)+ ...
                    (((y(:,1) + v_0_y) - repmat_Ym_0).^2)+ ...
                    (((z(:,1) + v_0_z) - repmat_Zm_0).^2)- ...
                    (R_0^2);

        % Vektor psi
        psi_0 = b_0;
        clear b_0;

        % Widerspruchsvektor
        w = -B * v_0 + psi_0;
        clear psi_0;

        % Normal equations
        N = [B*Q_ll*B',A;A',zeros(4,4)];

        rechte_seite = [-w;0;0;0;0];

        % Solving equation system
        delta_xk = N\rechte_seite;

        % Vektor delta_x extrahieren
        delta_x(1,1) = delta_xk((nx)+1,1);
        delta_x(2,1) = delta_xk((nx)+2,1);
        delta_x(3,1) = delta_xk((nx)+3,1);
        delta_x(4,1) = delta_xk((nx)+4,1);

        %Korrelatenvektor_k extrahieren
        k = zeros(nx,1);
        k(:,1) = delta_xk(1:nx,1);

        % Computation of unknown values
        Xm_dach = Xm_0 + delta_x(1,1);
        Ym_dach = Ym_0 + delta_x(2,1);
        Zm_dach = Zm_0 + delta_x(3,1);
        R_dach  = R_0  + delta_x(4,1);

        % Computation of residuals
        v_dach = Q_ll*B'*k;

        % Computed unknowns as new initial values
        Xm_0 = Xm_dach;
        Ym_0 = Ym_dach;
        Zm_0 = Zm_dach;
        R_0  = R_dach;

        % Computed residuals as new initial values
        v_0   = v_dach;
        v_0_x = v_dach(1       : nx  , 1);
        v_0_y = v_dach(nx   +1 : 2*nx, 1);
        v_0_z = v_dach(2*nx +1 : 3*nx, 1);

        max_delta_x = max(abs(delta_x));

        iteration = iteration + 1;

        if iteration > maxIterations
            disp('No convergation of LS');
            break;
        end

    end% --> Ende der Iterationsschleife

    % Computation of variance factor
    P = inv(Q_ll);
    vTPv      = v_dach'* P * v_dach;
    redundanz = n - length(delta_xk);

    sigma_0_a_post = sqrt(vTPv/redundanz);
        
    % Standard deviation of unknowns
    Ninv = inv(N);
    Q_xx = -Ninv(nx+1 : nx+4, nx+1 : nx+4);

    m_Xm  = sigma_0_a_post*(sqrt(Q_xx(1,1)));
    m_Ym  = sigma_0_a_post*(sqrt(Q_xx(2,2)));
    m_Zm  = sigma_0_a_post*(sqrt(Q_xx(3,3)));
    m_R   = sigma_0_a_post*(sqrt(Q_xx(4,4)));

    Xm_dach = Xm_dach + gravityC(1);
    Ym_dach = Ym_dach + gravityC(2);
    Zm_dach = Zm_dach + gravityC(3);

    % NV
    Q_kk_dach = Ninv(1:nx,1:nx);
    Q_vv = Q_ll * B' * Q_kk_dach * B * Q_ll;
    NV = abs(v_dach) ./ (sigma_0_a_pri*sqrt(diag(Q_vv)));
    
      %Least squares testing  
    
    %chi-square test
%     chiSquare = vartest(quality.sigma0_apost^2,sigma0_apri^2)
    f_quantil    = finv(confidenceLevel,redundanz, redundanz);
    f_testValue  = sigma_0_a_post^2/(redundanz * sigma_0_a_pri^2);
    fTest_result = f_testValue < f_quantil
    
     
    if dataSnooping == 1
    else
  
        iSnoop
        %t-test a post
        %   Unkorrelierte Beobachtung
        t_quantil = tinv(confidenceLevel,redundanz-1);
        sigma_0_a_post_strich = sqrt((vTPv - (v_dach .* v_dach)./diag(Q_vv))'/redundanz-1);
        t_testValue  = v_dach ./ (sigma_0_a_post_strich'.*sqrt(diag(Q_vv)));
        groberFehler = t_testValue > t_quantil;
        
%         groberFehler = find(NV > 3);

        if sum(groberFehler) == 0
            break;
        else   
            [sortGroberFehler sortIndexGroberFehler] = sort(NV(groberFehler),'descend');  
            maxiD = length(groberFehler);
            indexCheckedPoints = [];
            iD = 1;
            flagCheckPoints = true;
            while flagCheckPoints
                if groberFehler(sortIndexGroberFehler(iD)) > nx && groberFehler(sortIndexGroberFehler(iD)) <= nx * 2
                    indexCheckedPoints(end +1) =  groberFehler(sortIndexGroberFehler(iD)) -nx;
                elseif groberFehler(sortIndexGroberFehler(iD)) > (nx * 2)
                    indexCheckedPoints(end +1) =  groberFehler(sortIndexGroberFehler(iD)) -(nx * 2);
                else
                    indexCheckedPoints(end +1) =  groberFehler(sortIndexGroberFehler(iD));
                end

                indexCheckedPoints = unique(indexCheckedPoints);

                iD = iD +1;
                if iD > maxiD
                    flagCheckPoints = false;
                elseif length(indexCheckedPoints) >= dataSnooping
                    flagCheckPoints = false;
                end  
            end
            point3D(indexCheckedPoints,:) = [];
            Q_ll([indexCheckedPoints (indexCheckedPoints +nx) (indexCheckedPoints +nx*2)],:) = []; 
            Q_ll(:,[indexCheckedPoints (indexCheckedPoints +nx) (indexCheckedPoints +nx*2)]) = []; 
        end
    end
    iSnoop = iSnoop +1;
    
end % -> end datasnooping


kugel = struct('center',[Xm_dach Ym_dach Zm_dach],...
               'sigma_center', [m_Xm m_Ym m_Zm]  ,...
               'radius', R_dach                  ,...
               'sigma_radius', m_R               ,...
               'residuals',v_dach                ,...
               'Qxx',Q_xx                        ,...
               'sigma0_apost',sigma_0_a_post     ,...
               'iterations', iteration           ,...
               'balancedCoord', [x + v_0_x + gravityC(1), y + v_0_y + gravityC(2), z + v_0_z + gravityC(3)]);

%% Plot
if plot
   
    res_v_dach = sqrt(v_dach(1:nx) .* v_dach(1:nx) + v_dach(nx +1:nx*2) .* v_dach(nx +1:nx*2) + v_dach(2*nx +1:nx*3) .* v_dach(2*nx +1:nx*3));
    figure_res = figure();
    scatter3(x+gravityC(1),y+gravityC(2),z+gravityC(3),3,res_v_dach,'filled');
    colormap(jet(100)); 
    colorbar('location','EastOutside');
    axis equal;
    title ('Residuals resulting');
    xlabel ('X [m]');
    ylabel ('Y [m]');
    zlabel ('Z [m]');
    hold on;
    n = round(R_dach*100);
    [xi,yi,zi] = ellipsoid(Xm_dach,Ym_dach,Zm_dach,R_dach,R_dach,R_dach,n);
    [xi,yi,zi] = ellipsoid(Xm_dach,Ym_dach,Zm_dach,R_dach,R_dach,R_dach);
    h=surf(xi,yi,zi,'faceAlpha', 0.5,'EdgeColor','none');
    
    set(h,'FaceLighting','phong','FaceColor','interp',...
          'AmbientStrength',2)
    light('Position',[1 0 0],'Style','infinite');
    hold off;
    

    figure_vX = figure();
    title 'Residual resulting'
    scatter3(x+gravityC(1),y+gravityC(2),z+gravityC(3),3,v_dach(1:nx),'filled');
    colormap(jet(100)); 
    colorbar('location','EastOutside')
    axis equal
    title ('Residuals X');
    xlabel ('X [m]');
    ylabel ('Y [m]');
    zlabel ('Z [m]');

    figure_vY = figure();
    title 'Residual resulting'
    scatter3(x+gravityC(1),y+gravityC(2),z+gravityC(3),3,v_dach(nx +1:nx*2),'filled')
    colormap(jet(100)); 
    colorbar('location','EastOutside')
    axis equal;
    title ('Residuals Y');
    xlabel ('X [m]');
    ylabel ('Y [m]');
    zlabel ('Z [m]');

    figure_vZ = figure();
    title 'Residual resulting'
    scatter3(x+gravityC(1),y+gravityC(2),z+gravityC(3),3,v_dach(2*nx +1:nx*3),'filled')
    colormap(jet(100)); 
    colorbar('location','EastOutside')
    axis equal;
    title ('Residuals Z');
    xlabel ('X [m]');
    ylabel ('Y [m]');
    zlabel ('Z [m]');
    
end
%% Protokoll

% Final results
if exist('reportPath', 'var') 
    
    erg_ptr = fopen(reportPath,'w');

    fprintf(erg_ptr,'\n\t\t2D KUGELAUSGLEICHUNG 3D ERGEBNISPROTOKOLL');
    fprintf(erg_ptr,'\n\t\t-------------------------------------');
    fprintf(erg_ptr,'\n\t\t (Ausgleichung im GH-Modell, streng)');
    fprintf(erg_ptr,'\n\n\tBearbeiter    : Claudius Schmitt');

    ctime = clock;
    fprintf(erg_ptr,'\n\tBerechnet am  : %s um %d:%d Uhr',date,ctime(4),ctime(5));

    fprintf(erg_ptr,'\n\n\tsigma_0 a priori    :%9.4f',sigma_0_a_pri);
    fprintf(erg_ptr,'\n\tsigma_0 a posteriori:%9.4f',sigma_0_a_post);

    fprintf(erg_ptr,'\n\n\n\tAnzahl Beobachtungen: %7d',n);

    fprintf(erg_ptr,'\n\n\n\tVerbesserungen');
    fprintf(erg_ptr,'\n\t---------------\n');
    fprintf(erg_ptr,'\n\t\t\t      vx          vy          vz');
    fprintf(erg_ptr,'\n\t--------------------------------------------------');


    fprintf(erg_ptr,'\n\tmax:\t%10.4f  %10.4f  %10.4f  ',max(v_dach(1:nx,1)),max(v_dach(nx:2*nx,1)),max(v_dach(2*nx:3*nx,1)));
    fprintf(erg_ptr,'\n\tmin:\t%10.4f  %10.4f  %10.4f  ',min(v_dach(1:nx,1)),min(v_dach(nx:2*nx,1)),min(v_dach(2*nx:3*nx,1)));

    fprintf(erg_ptr,'\n\n\n\tUnbekannte');
    fprintf(erg_ptr,'\n\t----------');
    fprintf(erg_ptr,'\n\n\tKreisparameter');
    fprintf(erg_ptr,'\n\tXm: %12.6f  m.F.:%7.6f',Xm_dach,m_Xm);
    fprintf(erg_ptr,'\n\tYm: %12.6f  m.F.:%7.6f',Ym_dach,m_Ym);
    fprintf(erg_ptr,'\n\tZm: %12.6f  m.F.:%7.6f',Zm_dach,m_Zm);
    fprintf(erg_ptr,'\n\tR : %12.6f  m.F.:%7.6f',R_dach,m_R);

    fprintf(erg_ptr,'\n\n\n\tAnzahl der Iterationen (mit epsilon= %.4e): %d',EPSILON,iteration);
    if iteration == maxIterations
        fprintf(erg_ptr,'\n\n\n\tNo convergent in the given bounds of iteration');
    end;
    if dataSnooping == 0
        fprintf(erg_ptr,'\n\n\n\tNo DataSnooping');
    else
        fprintf(erg_ptr,'\n\n\n\tDataSnooping iterations: %2d',iSnoop);
    end
    fclose(erg_ptr);
    
else
end

%% time

ElapsedTimeSphereApprox = toc(tic_sphere)
