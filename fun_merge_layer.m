function [B1_pos_star,B2_neg_star] = fun_merge_layer(B1_pos,B1_neg,B2_pos, B2_neg, Phi, Lambda)
% fun_merge_layer  Merge the upmost layer into the boundary condition
%                  matrix.
%
% [B1_pos_star,B2_neg_star] = fun_merge_layer(B1_pos,B1_neg,B2_pos, B2_neg, Phi, Lambda)
%
% Input Parameters:
%   B1_pos,B1_neg: Boundary condition matrices at upper interface
%   B2_pos,B2_neg: Boundary condition matrices at lower interface
%   Phi, Lambda:   Decomposed tranfer matrix - [T] = [Phi][Lambda][Phi]^-1
% Output Parameters:  
%   B1_pos_star:   New BC matrix at '+' side of interface
%   B1_neg_star:   New BC matrix at '-' side of interface
%
% Ref: 
% [1] Song, Guochenhao, Zhuang Mo, and J. Stuart Bolton, "A Transfer-Matrix
%     -Based Approach to Predicting Acoustic Properties of a Layered System 
%     in a General, Efficient, and Stable Way," SAE Int. J. Adv. & Curr. 
%     Prac. in Mobility 6(2):922-934, 2024.
%
% Written by: 
% Guochenhao Song
% Ray W. Herrick Lab, Purdue University
% Email: guochenhaosong@gmail.com
% 2022 Fall

% dimension for new matrices
N = size(Phi,1);
N1 = size(B1_pos,1)+size(B2_neg,1);
N2 = size(B1_pos,2)+size(B2_neg,2)+N;
row_st = size(B1_pos,1)+1;
col_st = size(B1_pos,2)+1;
col_ed = size(B1_pos,2)+N;

% construct big matrix
A = zeros(N1, N2);
A(1:row_st-1, 1:col_st-1)    = B1_pos;
A(1:row_st-1, col_st:col_ed) = -B1_neg*Phi*Lambda;
A(row_st:end, col_st:col_ed) = B2_pos*Phi;
A(row_st:end, col_ed+1:end)  = -B2_neg;

% Elimination
% A_BC = A(N+1:end,col_st:col_ed)/A(1:N,col_st:col_ed)*A(1:N,:);
A_BC = A(N+1:end,col_st:col_ed)*inv(A(1:N,col_st:col_ed))*A(1:N,:);

% extract new boundary condition matrices
B1_pos_star = A(N+1:end,1:col_st-1) - A_BC(:, 1:col_st-1);
B2_neg_star = -(A(N+1:end,col_ed+1:end) - A_BC(:, col_ed+1:end));
end


