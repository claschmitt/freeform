%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Basis function
%
% compute parameters for global surface interpolation
% INPUT
% r = control point index function with uk starting with 1 - columns
% s = control point index function with vl starting with 1 - rows
% Q = control Points as a 2 dimensional Cell array of a grided point net
%     where each Cell consist of an 3 dim array [x y z]
% OUTPUT
% uk = value in the horizontal knot vector
% vl = value in the vertical knot vector
% 
% cs, 22.05.2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [uk, vl] = SurfMeshParams (r,s,Q)

% compute uk
% ------------------
numU = s + 1; %number of rows
uk = zeros(1,r +1);
uk(0 +1) = 0.0;
uk(r +1) = 1.0;

for l = 0 :  s
    total = 0.0;
    cds = zeros(r +1,1);
        
    for k = 1 : r
        index_p1 = s*l+k-1 +1;
        index_p2 = s*l+k +1;
        cds(k +1) = distPoint2Point(Q(index_p1,:), Q(index_p2,:));
        total = total + cds(k +1);
    end
        
    if total == 0.0
            numU = numU - 1;
    else
        %Eq 9.8
        d = 0.0;
        for k = 1 : r-1
            d = d + cds (k +1);
            uk(k +1) = uk(k +1) + d/total;
        end
    end

    if numU == 0
        msgbox ('Surfmesh error');
    end
end
    
uk = uk/numU;
uk(1,r +1) = 1.0;

% ----------------------------------------------------------------------
% compute vl
if s == 0 
    vl = 0;
else
    numV = r + 1; %number of rows
    vl = zeros(1,s +1);
    vl(0 +1) = 0.0;
    vl(s +1) = 1.0;



    for k = 0 :  r
        total = 0.0;
        cds = zeros(s +1,1);

        for l = 1 : s
            index_p1 = k+(r+1)*(l-1) +1;
            index_p2 = k+r+1+(r+1)*(l-1) +1;
            cds(l +1) = distPoint2Point(Q(index_p1,:), Q(index_p2,:));
            total = total + cds(l +1);
        end

        if total == 0.0
                numV = numV - 1;
        else
            %Eq 9.8
            d = 0.0;
            for l = 1 : s-1
                d = d + cds (l +1);
                vl(l +1) = vl(l +1) + d/total;
            end
        end

        if numV == 0
            msgbox ('Surfmesh error');
        end
    end

    vl = vl/numV;
    vl(s +1) = 1.0;
end

