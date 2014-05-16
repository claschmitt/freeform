%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot curve
%
% INPUT
% nurbs = nurbs struct
% subInterval = subsample rate of the curve
% figureCurve = figure objekt
% axesHandle = axes Handle
% OUTPUT
% plotting nurbs curve
%
% cs, 22.08.2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function subPoints = plotCurve(nurbs,subInterval,figureCurve)
% subInterval=100;

[r, dimPoints] = size(nurbs.coefs);
subPosition = linspace(0.0,1.0,subInterval);
subPoints = zeros(subInterval,dimPoints);
figure(figureCurve);

for i=1:subInterval
    subPoints(i,:) = curvePoint(subPosition(i), nurbs);
end

% % Plot curve
% plot3(subPoints(:,1),subPoints(:,2),subPoints(:,3), '-b');
% hold on;
% % Plot controle polygon
% plot3(nurbs.coefs(1,:),nurbs.coefs(2,:), nurbs.coefs(3,:),'--r');

switch dimPoints
    case 2
        plot(subPoints(:,1),subPoints(:,2), '-b');
        hold on;
        % Plot controle polygon
        plot(nurbs.coefs(:,1), nurbs.coefs(:,2),'--.r');
        
        hold off;
    case 3
        plot3(subPoints(:,1), subPoints(:,2),subPoints(:,3), '-b');
        hold on;
        % Plot controle polygon
        plot3(nurbs.coefs(:,1), nurbs.coefs(:,2), nurbs.coefs(:,3),'--.r');
        
        hold off;
    otherwise
        disp('ERROR - points Dimesion')
end

