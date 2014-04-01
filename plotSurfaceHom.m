%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot Surface
%
% plot the aproximated nurbs surface
% input
% nurbs = struct of nurbsparameters
% subInterval = sampling rate of plottet points
% surfacePolygon = Boundingbox of Surface
% measPoint = measured points
% meanPoint = 
% output
% [X Y Z] = plotted point
% 
% cs, 12.02.2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [intPoint, intPoint_normal] = plotSurfaceHom(nurbs,subInterval, figureSurface, surfacePolygon, measPoint, quality)

%% testing 


%%
subPosition = linspace(0.0,1.0,subInterval);

intPoint = zeros(subInterval, subInterval,3);
intPoint_row = zeros (subInterval,3);
dersS_newX  = zeros (subInterval,subInterval);
dersS_newY  = zeros (subInterval,subInterval);
dersS_newZ  = zeros (subInterval,subInterval);

for i=1:subInterval
    b = 1;
    tmp_intPoint  = zeros(subInterval,3);
    dersS     = zeros(subInterval,3);
    
    for j=1 : subInterval

        [tmp_intPoint] = surfacePoint(nurbs.numberU -1, nurbs.orderU -1, nurbs.knotsU,nurbs.numberV -1, nurbs.orderV -1, nurbs.knotsV, nurbs.coefs',subPosition(i),subPosition(j));
       
        if isempty(tmp_intPoint)
        else
            intPoint_row(b,:) = tmp_intPoint(1,:);
            b = b+1;
        end
    end
    intPoint (:,i,1) = intPoint_row(:,1);
    intPoint (:,i,2) = intPoint_row(:,2);
    intPoint (:,i,3) = intPoint_row(:,3);

    
    dersS_newX (:,i) = dersS(:,1);
    dersS_newY (:,i) = dersS(:,2);
    dersS_newZ (:,i) = dersS(:,3);
    
end


%% plotting

figure(figureSurface);

% 	Select plotting points
subplot(2,2,1);

colorbarScale = [0.0200 0.0150 0.0100 0.0075 0.0050 0.0025 0.0000 -0.0025 -0.0050 -0.0075 -0.0100 -0.0150 -0.0200];

%plot surface with residuals
subplot(3,2,1);
res_quality = sqrt(quality.residuals(:,1) .* quality.residuals(:,1) +  quality.residuals(:,2) .* quality.residuals(:,2) + quality.residuals(:,3) .* quality.residuals(:,3));
resSurface = griddata(measPoint(:,2),measPoint(:,1),measPoint(:,3),res_quality,intPoint(:,:,2),intPoint(:,:,1),intPoint(:,:,3), 'nearest');
surf(intPoint(:,:,2),intPoint(:,:,1),intPoint(:,:,3),resSurface, 'EdgeColor','none','LineStyle','none','FaceLighting','phong');
colormap(jet(100)); 
caxis([ min(colorbarScale)  max(colorbarScale)]);
colorbar('location','EastOutside','yTick', colorbarScale(1:7), 'YTickLabel', colorbarScale(1:7));
shading faceted;
axis equal;
title ('Surface with resulting residuals');
grid on;
xlabel ('X [m]');
ylabel ('Y [m]');
zlabel ('Z [m]');

subplot(3,2,2);
surf(intPoint(:,:,2),intPoint(:,:,1),intPoint(:,:,3)); 
shading faceted;
axis equal;
title ('Surface');
grid on;
xlabel ('X [m]');
ylabel ('Y [m]');
zlabel ('Z [m]');
hold on;


%plot normal vektors of each point
[intPoint_normal(:,:,1) intPoint_normal(:,:,2) intPoint_normal(:,:,3)] = surfnorm(intPoint(:,:,1), intPoint(:,:,2), intPoint(:,:,3));

quiver3(intPoint(1:10:end,1:10:end,2),intPoint(1:10:end,1:10:end,1),intPoint(1:10:end,1:10:end,3), ...
    intPoint_normal(1:10:end,1:10:end,2), intPoint_normal(1:10:end,1:10:end,1), intPoint_normal(1:10:end,1:10:end,3),1, 'color' , 'blue');
hold off


%plot measured points with residuals
numPlotingPoints = 10000;
if length(measPoint) > numPlotingPoints
    randIndex_measPoint = randi(length(measPoint),numPlotingPoints,1);
else
    randIndex_measPoint = ones(length(measPoint),1);
end

subplot(3,2,3);
scatter3(measPoint(randIndex_measPoint,2), measPoint(randIndex_measPoint,1), measPoint(randIndex_measPoint,3), 2.2, res_quality(randIndex_measPoint), 'filled');
colormap(jet(100)); 
caxis([ min(colorbarScale)  max(colorbarScale)]);
colorbar('location','EastOutside','yTick', colorbarScale(1:7), 'YTickLabel', colorbarScale(1:7));
axis equal;
title ('Residuals resulting');
xlabel ('X [m]');
ylabel ('Y [m]');
zlabel ('Z [m]');

subplot(3,2,4);
title 'Residual resulting'
scatter3(measPoint(randIndex_measPoint,2), measPoint(randIndex_measPoint,1), measPoint(randIndex_measPoint,3), 2.2, quality.residuals(randIndex_measPoint,1), 'filled');
colormap(jet(100)); 
caxis([min(colorbarScale)  max(colorbarScale)]);
colorbar('location','EastOutside','yTick', colorbarScale, 'YTickLabel', colorbarScale);
axis equal
title ('Residuals X');
xlabel ('X [m]');
ylabel ('Y [m]');
zlabel ('Z [m]');

subplot(3,2,5);
title 'Residual resulting'
scatter3(measPoint(randIndex_measPoint,2), measPoint(randIndex_measPoint,1), measPoint(randIndex_measPoint,3), 2.2, quality.residuals(randIndex_measPoint,2), 'filled');
colormap(jet(100)); 
caxis([min(colorbarScale)  max(colorbarScale)]);
colorbar('location','EastOutside','yTick', colorbarScale, 'YTickLabel', colorbarScale);
axis equal;
title ('Residuals Y');
xlabel ('X [m]');
ylabel ('Y [m]');
zlabel ('Z [m]');

subplot(3,2,6);
title 'Residual resulting'
scatter3(measPoint(randIndex_measPoint,2), measPoint(randIndex_measPoint,1), measPoint(randIndex_measPoint,3), 2.2, quality.residuals(randIndex_measPoint,3), 'filled');
colormap(jet(100)); 
caxis([min(colorbarScale)  max(colorbarScale)]);
colorbar('location','EastOutside','yTick', colorbarScale, 'YTickLabel', colorbarScale);
axis equal;
title ('Residuals Z');
xlabel ('X [m]');
ylabel ('Y [m]');
zlabel ('Z [m]');

hold off;
