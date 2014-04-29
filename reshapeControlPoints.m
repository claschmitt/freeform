%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reshapeControlPoints
%
% Resphape control point matrix to point parameters
% INPUT
%   pointsParameter = [uk vl];
%   nurbs  = struct('form','', 'dimU','', 'dimV','', 'numberU','', 'numberV','', 'coefs','', 'orderU','', 'orderV','', 'knotsU','', 'knotsV','');
%  optional parameters
%   quality = quality structur of Apporximation especially with the Qxx
%
% OUTPUT
%   CPMatrix = Control point Matrix to pointsParameter for Equation
%  optional
%   Q_CPMatrix = Kovariance Matrix of the control points for Equation
%
%
% cs, 25.04.2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [CPMatrix, varargout] = reshapeControlPoints(pointsParameter, nurbs, varargin)

if isempty(varargin)
    varargout = [];
    flagQuality = false;
else
    quality = varargin{1};
    flagQuality = true;
    Qss = zeros(nurbs.orderU * nurbs.orderV * 3);
end

spanU = findSpan(pointsParameter(1),nurbs.knotsU);
spanV = findSpan(pointsParameter(2),nurbs.knotsV);

CPMatrix = zeros(nurbs.orderU * nurbs.orderV * 3,1);

for iterBasisU = 0 : nurbs.dimU
    
    % filling CP Matrix with basis functions
    indexCoefs_column_start = int64((((iterBasisU + spanU - nurbs.orderU) * nurbs.numberV)) +  (spanV - nurbs.orderV) +1);
    indexCoefs_column_end   = int64((((iterBasisU + spanU - nurbs.orderU) * nurbs.numberV)) +  (spanV - nurbs.orderV) + nurbs.orderV);
    
    indexCPMatrix_row_start = iterBasisU * nurbs.orderV +1;
    indexCPMatrix_row_end   = iterBasisU * nurbs.orderV + nurbs.orderV; 
    
    CPMatrixX(indexCPMatrix_row_start : indexCPMatrix_row_end, 1) = nurbs.coefs(indexCoefs_column_start : indexCoefs_column_end,1);
    CPMatrixY(indexCPMatrix_row_start : indexCPMatrix_row_end, 1) = nurbs.coefs(indexCoefs_column_start : indexCoefs_column_end,2);
    CPMatrixZ(indexCPMatrix_row_start : indexCPMatrix_row_end, 1) = nurbs.coefs(indexCoefs_column_start : indexCoefs_column_end,3);
    
    if flagQuality
           
        indexQxx_start (iterBasisU +1) = int64((((iterBasisU + spanU - nurbs.orderU) * nurbs.numberV * 3)) +  (spanV - nurbs.orderV) +1);
        indexQxx_end   (iterBasisU +1) = int64((((iterBasisU + spanU - nurbs.orderU) * nurbs.numberV * 3)) +  (spanV - nurbs.orderV) + (nurbs.orderV * 3));
        
%         Qss(iterBasisU * nurbs.orderU * 3 +1, :) = quality.Qxx(indexQss_row_start : indexQss_row_end, indexQss_column_start : indexQss_column_end)
%         Qss(iterBasisU * nurbs.orderU * 3 +2, :) = quality.Qxx(indexQss_row_start : indexQss_row_end, indexQss_column_start : indexQss_column_end)
%         Qss(iterBasisU * nurbs.orderU * 3 +3, :) = quality.Qxx(indexQss_row_start : indexQss_row_end, indexQss_column_start : indexQss_column_end)
    end
end

% Qss = [ quality.Qxx(indexQxx_start(1,1) : indexQxx_end (1,1), indexQxx_start(1,1) : indexQxx_end (1,1)) quality.Qxx(indexQxx_start : indexQxx_end, indexQxx_start : indexQxx_end)


CPMatrix(1:3:end,1) = CPMatrixX;
CPMatrix(2:3:end,1) = CPMatrixY;
CPMatrix(3:3:end,1) = CPMatrixZ;

if flagQuality
    varargout{1} = Qss;
else
    varargout = [];
end