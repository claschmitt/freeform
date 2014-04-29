%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Global Curve Approximation
%
% compute least squares fitting of a non uniform B-Spline curve with fixed control points
% input
% r = max point index of measured points in each row
% Q = mesured points
% p = degree of function 1 with Nots U
% n = max index of number of control Points (U-direction)starting with 0
% U = knot vektor of function 1 with degree p
% P = Control Points of Controle net
% 
% cs, 25.04.2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [U,P_full] = globalCurveApprox (r,Q,p,n,uk,U,Q_ll)


% compute internal knot vector U
if U<=0
    du = (r + 1) / (n - p + 1);

    U = zeros((p+1)*2+n-p,1);
    U (p+1+n-p +1:(p+1)*2+n-p)= 1; 

    for k=1 : n-p
       i = floor(k * du);
       alpha = k * du - i;
       U(p+k +1) = (1 - alpha) * uk(i-1 +1) + alpha * uk(i +1);
    end
end


% fill coefficient matrix in U direction without the first and last points 
Au_full = zeros((r +1), (n +1));
pointsRu = zeros(n-1,3);
indexRu_column = 1;
RuX = zeros(r-1,1);
RuY = zeros(r-1,1);
RuZ = zeros(r-1,1);

for k=0 : r 
%         span = findspan(n+p+1,p,uk(k+1 +1),U);
% start and end function compliant with span index
    if uk(k +1) >= U(n+1 +1)
        span = length(U)-p-1 -1;
    else
        tmp_span = find(U<=uk(k +1));
        span = max(tmp_span)-1;
    end

%     N = basisFunction(span,U(span +1),p,U);
    N = basisFunction(span,uk(k +1),p,U);

    if n<p
        msgbox('number of controle points must be higher','error');
        break;
    end 


    indexAu_column_start = int64(span-p +1);
    indexAu_column_end   = int64(span-p+p +1);
    indexN_start         = int64(0 +1);
    indexN_end           = int64(p +1);  %not clear yet


    indexAu_row = k +1;
    Au_full(indexAu_row,indexAu_column_start : indexAu_column_end) = N(indexN_start : indexN_end);

    if k>=1 && k<=r 
        index_Q0 = 0 +1;
        index_Qr = r +1;
%         pointsRu(k,:) = Q(k +1,:) -  Au_full(k +1,0 +1).*Q(index_Q0,:) - Au_full(k +1,n +1).*Q(index_Qr,:);
        %         pointsRu(k,:) = Q(k +1,:);
    end

end



% 
% Au(:,:) = Au_full((1 +1):(r-1 +1), (1 +1):(n-1 +1));

% Ru = zeros((n-1),3);

% filling Ru Matrix
% for i=1 : n-1
%     Ru(i,1) = sum(Au(:,i) .* pointsRu(:,1));
%     Ru(i,2) = sum(Au(:,i) .* pointsRu(:,2));
%     Ru(i,3) = sum(Au(:,i) .* pointsRu(:,3));
% end

% Au = sparse(Au);
% Ru = sparse(Ru);

% P(:,1) = (Au' / Q_ll(3:2:end-2, 3:2:end-2) * Au) \ Ru(:,1);
% P(:,2) = (Au' * Au) \ Ru(:,2);
% P(:,3) = (Au' / Q_ll(4:2:end-2,  4:2:end-2) * Au) \ Ru(:,3);

% % Reshape Matrixes
Q_ll = sparse(Q_ll);
Au_full = sparse(Au_full);
Au_xz = zeros(length(Q_ll),2*(n+1));
Au_xz(1:2:end, 1:2:end) = Au_full;
Au_xz(2:2:end, 2:2:end) = Au_full;

pointsRuTMP = reshape([Q(:,1) Q(:,3)]',length(Q) * 2,1);

P = (Au_xz' / Q_ll * Au_xz) \ (Au_xz' / Q_ll * pointsRuTMP);
pointsRuTMP_dach = Au_xz * P;
v_dach = pointsRuTMP_dach - pointsRuTMP;
varianz_dach = (v_dach' * v_dach)/(length(pointsRuTMP)-n+1);
sigma0_dach = sqrt(varianz_dach)


P_full = zeros(n+1,3);
% P_full(1,:) = Q(1,:);
% P_full(n+1,:) = Q(r+1,:);
P_full(:,1) = P(1:2:end,1); 
% P_full(2:n,2) = P(:,);
P_full(:,3) = P(2:2:end,1);

% plot(Q(1:r+1,1), Q(1:r+1,3),'o', 'color', 'blue');
% hold on;
% plot(P_full(:,1),P_full(:,3),'-o','color', 'red')
% hold off;
