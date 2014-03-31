%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IncidenceAngle
%
% Calculate the Incidence angle of one direction within the object
% structure
%
% input
% stationPoint = Station coordinates [y z]
% MeasuredPoints [y z]
% geometryFlag
%   1 = straight line

% out
% iA = for each point in verticle direction 
% 
% cs, 10.02.2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function iA = IncidenceAngle(stationPoint, measuredPoints, geometryFlag)

%% test
% load('Y:\Forschung\AusIng\institut_fuer_massivbau\rethen\20130701_Messung 3\Scans\Matlab\Session2\Nr14183 ultrahigh I637 Profiler.part10.xyz.asc.mat');
% measuredPoints = scan.data{1}(:,2:3);
% measuredPoints(:,1) = measuredPoints(:,1)  - abs(min(scan.boundingBox(:,1)));
% clear scan;
% stationPoint = [0 0];
% geometryFlag = 1;

%%
switch geometryFlag
    case 1
        [p,S] = polyfit(measuredPoints(:,1), measuredPoints(:,2),1);
        
        observations = zeros(length(measuredPoints) * 2,1);
        observations(1:2:end,1) = measuredPoints(:,1);
        observations(2:2:end,1) = measuredPoints(:,2);
        clear measuredPoints;
    
        % rotate koordinate system to y along straingth line
        rotationAngle = atan(p(1));
        Rx = [  cos(rotationAngle) sin(rotationAngle) ;...
               -sin(rotationAngle) cos(rotationAngle) ];
           
        mainDiag(1:length(observations),1) = Rx(1,1);
%         mainDiag = diag(mainDiag);
        R = sparse(diag(mainDiag));
        clear mainDiag;
        
        minorDiagAbove = zeros(length(observations)-1,1);
        minorDiagAbove (1:2:end) = Rx(1,2);
%         minorDiagAbove = diag(minorDiagAbove,1);
        R = sparse(R + diag(minorDiagAbove,1));
        clear minorDiagAbove;
        
        minorDiagBelow = zeros(length(observations)-1,1);
        minorDiagBelow (1:2:end) = Rx(2,1);
%         minorDiagBelow = diag(minorDiagBelow,-1);
        R = sparse(R + diag(minorDiagBelow,-1));
        clear minorDiagBelow;
        
        observationsRot = R * observations;
        clear observations;
        stationPointRot = stationPoint * Rx;
        clear stationPoint;
       
        % process direction angle
        riWi = atan2(observationsRot(1:2:end,1) - stationPointRot(1) , observationsRot(2:2:end,1) - stationPointRot(2));
        iA = riWi;
        
        % Quadrant query 
        % find(riWi < (pi/2))
%         indexQ2 = find(riWi > (pi/2));
%         indexQ3 = find(riWi > (pi  ));
%         indexQ4 = find(riWi > (pi + pi/2));
% 
%         iA(indexQ2) = riWi(indexQ2) - pi/2;
%         iA(indexQ3) = riWi(indexQ3) - pi;
%         iA(indexQ4) = riWi(indexQ4) - pi - pi/2;

end