%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute internal knot vector
%
% compute the internal knot vector to create the Berstein polynomials
% INPUT
% pointsParameter = mesured points in parameter form [uk] or [uk vl]
% nurbs           = nurbs structure
% method = 'uniform' or 'piegl_tiller'( for one or two dim, when
% methodV is not choosen)
% methodV = for second dimension calc method
%
% OUTPUT
% nurbs = nurbs filled with internal knot vector for the Berstein polynomial
%
% The NURBS book page 365 ff and 412
% cs, 02.05.2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function nurbs = internalKnots (pointsParameter,nurbs,method, methodV)


[numRow numCol] = size(pointsParameter);

if numRow == 1 && numCol > 1
    dim = 1;
    numPoints = numCol;
    pointsParameters = pointsParameters';
elseif numRow == 2 && numCol > 2
    dim = 2;
    numPoints = numCol;
    pointsParameters = pointsParameters';
elseif numCol == 1 && numRow > 1
    dim = 1;
    numPoints = numRow;
elseif numCol == 2 && numRow > 2
    dim = 2;
    numPoints = numRow;
else
    disp('Error in pointsParameter');
    quit;
end

if nargin == 3 && dim == 2
    methodV = method;
elseif dim == 1
    methodV = 0;
end

switch method
    case 'uniform'
        nurbs.knotsU = zeros(nurbs.orderU + nurbs.numberU,1);
        nurbs.knotsU(nurbs.orderU + nurbs.numberU - nurbs.orderU +1 : nurbs.orderU + nurbs.numberU) = 1;
        
        numSpans = nurbs.numberU - nurbs.orderU *2 +1;
        for iterKnots = 1 : nurbs.numberU - nurbs.orderU *2
            nurbs.knotsU (nurbs.orderU + iterKnots) = 1 /  numSpans * iterKnots
        end
        
    case 'piegl_tiller'
        
        
        d = numPoints / (nurbs.numberU - nurbs.orderU +1);
        uk_sorted = sort(pointsParameter(:,1));
        
        nurbs.knotsU = zeros(nurbs.orderU + nurbs.numberU,1);
        nurbs.knotsU(nurbs.orderU + nurbs.numberU - nurbs.orderU +1 : nurbs.orderU + nurbs.numberU) = 1;
        
        for iterKnots=1 : nurbs.numberU - nurbs.orderU
            i = floor(iterKnots * d);
            alpha = iterKnots * d - i;
            nurbs.knotsU(nurbs.orderU +iterKnots) = (1 - alpha) * uk_sorted(i-1 +1) + alpha * uk_sorted(i +1);
        end
        
        clear uk_sorted;
        
    otherwise
        disp ('Error in knot calculation');
        quit;
end

switch methodV
    case 'uniform'
        nurbs.knotsV = zeros(nurbs.orderV + nurbs.numberV,1);
        nurbs.knotsV(nurbs.orderV + nurbs.numberV - nurbs.orderV +1 : nurbs.orderV + nurbs.numberV) = 1;
        
        numSpans = nurbs.numberV - nurbs.orderV *2 +1;
        for iterKnots = 1 : nurbs.numberV - nurbs.orderV *2
            nurbs.knotsV (nurbs.orderV + iterKnots) = 1 /  numSpans * iterKnots
        end
        
    case 'piegl_tiller'
        
        d = numPoints / (nurbs.numberV - nurbs.orderV +1);
        vl_sorted = sort(pointsParameter(:,2));
        
        nurbs.knotsV = zeros(nurbs.orderV + nurbs.numberV,1);
        nurbs.knotsV(nurbs.orderV + nurbs.numberV - nurbs.orderV +1 : nurbs.orderV + nurbs.numberV) = 1;
        
        for iterKnots=1 : nurbs.numberV - nurbs.orderV
            i = floor(iterKnots * d);
            alpha = iterKnots * d - i;
            nurbs.knotsV(nurbs.orderV +iterKnots) = (1 - alpha) * vl_sorted(i-1 +1) + alpha * vl_sorted(i +1);
        end
        
        clear vl_sorted;
end