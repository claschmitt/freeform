%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FillAGM
%
% Filling A matrix with basis functions for Gau�-Markov model
% INPUT
%   pointsCartesian = [x y z];
%   pointsParameter = [uk vl];
%   nurbs  = struct('form','', 'dimU','', 'dimV','', 'numberU','', 'numberV','', 'coefs','', 'orderU','', 'orderV','', 'knotsU','', 'knotsV','');
%  optional parameters
%
% OUTPUT
%   AMatrix = Full A Matrix for system equation
%
%
% cs, 25.04.2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function AMatrix = fillAGM (pointsCartesian, pointsParameter, nurbs)
tic_fill_A = tic;

%% testing parameter

if nurbs.numberU < nurbs.orderU
    msgbox('number of controle points must be higher than degree of basis function in U direction','error');
end

if nurbs.numberV < nurbs.orderV
    msgbox('number of controle points must be higher than degree of basis function in V direction','error');
end

%% dimension test reordering points
if sum(pointsCartesian(:,1)) ~= 0
    flagX = 1;
else
    flagX = 0;
    pointsCartesian(:,1) = [];
end

if sum(pointsCartesian(:,2)) ~= 0
    flagY = 1;
else
    flagY = 0;
    pointsCartesian(:,2) = [];
end

if sum(pointsCartesian(:,3)) ~=0
    flagZ = 1;
else
    flagZ = 0;
    pointsCartesian(:,3) = [];
end

if sum(pointsParameter(:,1)) >= 0
    flagUk = 1;
else
    flagUk = 0;
    pointsParameter(:,1) = [];
end

if sum(pointsParameter(:,2)) >= 0
    flagVk = 1;
else
    flagVk = 0;
    pointsParameter(:,2) = [];
end


dimPoints = flagX + flagY + flagZ;


%% filling A
if dimPoints == 2
    
elseif dimPoints == 3
    
    [numPointsCartesian, c] = size(pointsCartesian);
    
    if numPointsCartesian == 1 % for single point calculation
        A_full   = zeros(numPointsCartesian * dimPoints, nurbs.orderU * nurbs.orderV * dimPoints);
        A_single = zeros(numPointsCartesian,nurbs.orderU * nurbs.orderV);
    else % for complet GM model
        A_full   = zeros(numPointsCartesian * dimPoints, nurbs.numberU * nurbs.numberV * dimPoints);
        A_single = zeros(numPointsCartesian, nurbs.numberU * nurbs.numberV);
    end
    
    %   filling A fore one direction
    for iterPoints=1 : numPointsCartesian
        
        % find span at internal knot vector
        spanU = findSpan(pointsParameter(iterPoints,1),nurbs.knotsU);
        spanV = findSpan(pointsParameter(iterPoints,2),nurbs.knotsV);
        
        % calc basis functions
        Nu = basisFunction(spanU -1,pointsParameter(iterPoints,1),nurbs.orderU -1,nurbs.knotsU);
        Nv = basisFunction(spanV -1,pointsParameter(iterPoints,2),nurbs.orderV -1,nurbs.knotsV);
        
        tensor_u_v = kron(Nu',Nv);
        
        for iterBasisU = 1 : nurbs.orderU
            
            % filling A with basis functions
            if numPointsCartesian == 1
                indexA_column_start = int64(((iterBasisU -1) * nurbs.orderV)  +1);
                indexA_column_end   = int64(((iterBasisU -1) * nurbs.orderV) + nurbs.orderV);
            else
                indexA_column_start = int64((((iterBasisU -1 + spanU - nurbs.orderU) * nurbs.numberV)) +  (spanV - nurbs.orderV) +1);
                indexA_column_end   = int64((((iterBasisU -1 + spanU - nurbs.orderU) * nurbs.numberV)) +  (spanV - nurbs.orderV) + nurbs.orderV);
            end
            A_single(iterPoints,indexA_column_start : indexA_column_end) = tensor_u_v(iterBasisU, :);
        end
        
    end
      
    %   filling A with all directions (X,Y,Z)
    % X
    A_full(1:3:end, 1:3:end) = A_single;
    % Y
    A_full(2:3:end, 2:3:end) = A_single;
    % Z
    A_full(3:3:end, 3:3:end) = A_single;
    
    AMatrix = A_full;
%     spy(AMatrix,'k');
    
    
else
    disp('Dimension error');
end


ElapsedTimeFill_A = toc(tic_fill_A);