%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% createPointsParametersCurve
%
% compute the internal knot vector to create the Berstein polynomials
% INPUT
% pointsKartesian = [X Y Z]
% method = 'centirpetal', 'chord', 'uniform' 
%
% OUTPUT
% pointsParameter = [uk vl] 
%
% The NURBS book page 365 ff
% cs, 02.05.2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pointsParameters = createPointsParametersCurve(pointsKartesian, method)

numPoints             = length(pointsKartesian);
pointsParameters      = zeros(numPoints ,1);

switch method
    case 'uniform'
        
        pointsParameters(2,:) = (2 : (numPoints -1)) / (numPoints  -2);
                
    case 'chord'
        d_single = zeros(numPoints -2,1);
        d_sum    = zeros(numPoints -2,1);
        
        for iterPoints = 1 : numPoints -2
            d_single (iterPoints) = distpoint2point(pointsKartesian(iterPoints,:), pointsKartesian(iterPoints +1,:));
            if iterPoints == 1
                d_sum    (iterPoints) =  d_single (iterPoints);
            else
                d_sum    (iterPoints) = d_sum(iterPoints-1) + d_single (iterPoints);
            end
        end
        
        pointsParameters(2:end) = d_sum ./ d_sum(end);
        
    case 'centripetal'
        d_single = zeros(numPoints -2,1);
        d_sum    = zeros(numPoints -2,1);
        
        for iterPoints = 1 : numPoints -1
            d_single (iterPoints) = sqrt(distPoint2Point(pointsKartesian(iterPoints,:), pointsKartesian(iterPoints +1,:)));
            if iterPoints == 1
                d_sum    (iterPoints) =  d_single (iterPoints);
            else
                d_sum    (iterPoints) = d_sum(iterPoints-1) + d_single (iterPoints);
            end
        end
        
        pointsParameters(2:end) = d_sum / d_sum(end);
        
    otherwise
        disp ('Error in knot calculation');
        quit;
end
