%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Basis function
%
% compute the nonvanishing Basis functions for B-Spline
% input
% i = knot span index (position in knot vector)!knot vector startindex = 1 !
% u = explicit knot value
% p = degree of the function
% U = knot vektor startindex = 1
% hallo alle leute
% hallo2
%hall4o
%hallo0
%hallo4
% cs, 11.05.2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [N] = basisFunction (i,u,p,U)

N = zeros(1,p +1);

N(0 +1) = 1.0;
% if u == 1
%     N(p +1) = 1;
% else
    for j=1 : p
    %     differences of knots
        indexU_left  = i+1-j +1;
        indexU_right = i+j +1;
        left (j +1) = u - U(indexU_left);
        right(j +1) = U(indexU_right) - u;
        saved = 0.0;

        for r= 0 : j-1 
            temp = N(r +1) / (right(r+1 +1) + left (j-r +1));
            N(r +1) = saved + right(r+1 +1) * temp;
            saved = left(j-r +1) * temp;
        end
        N(j +1) = saved;
    end
% end

% %Check weather calculation is right
%  if sum(N(1:p +1)) ~ 1.0
% else
%     N (:) = [];
% end


% _________________________________________________________________________
% Formula starting with index by 1
% 
% function [N] = BasisFunction (i_1,u,p_1,U)
% 
% N = zeros([1,p_1]);
% % N(0) = 1.0;
% N(1) = 1.0;
% 
% % for j=1 : p
% for j=2 : p_1
%     left (j) = u - U(i_1+1-j);
%     right(j) = U(i_1+j)-u;
%     saved = 0.0;
%     
%     for r= 1 : j-1 %problem hier trap zero lest 
% %         temp = N(r) / (right(r+1) + left (j-r));
% %         N(r) = saved + right(r+1) * temp;
% %         saved = left(j-r) * temp;
%         temp = N(r) / (right(r+1) + left (j-r+1));
%         N(r) = saved + right(r+1) * temp;
%         saved = left(j-r+1) * temp;
%     end
%     N(j) = saved;
% end
% 
% %Check weather calculation is right
% if sum(N(1:p_1)) ~ 1.0
% else
%     N = 'NaN';
% end
%     
