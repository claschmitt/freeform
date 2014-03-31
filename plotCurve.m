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

subPosition = linspace(0.0,1.0,subInterval);
subPoints = zeros(subInterval,3);
figure(figureCurve);
    
for i=1:subInterval
    subPoints(i,:) = curvePoint(nurbs.numberU -1, nurbs.orderU -1, nurbs.knotsU, nurbs.coefs',subPosition(i));
end

% % Plot curve
% plot3(subPoints(:,1),subPoints(:,2),subPoints(:,3), '-b');  
% hold on;
% % Plot controle polygon
% plot3(nurbs.coefs(1,:),nurbs.coefs(2,:), nurbs.coefs(3,:),'--r');

plot(subPoints(:,1),subPoints(:,3), '-b');  
hold on;
% Plot controle polygon
plot(nurbs.coefs(1,:), nurbs.coefs(3,:),'--.r');

hold off;

