function [TM] = fun_1layer_pred(B1_pos,B1_neg,B2_pos, B2_neg, Phi, Lambda)
% fun_1layer_pred  Solve for the 2x2 transfer matrix of a 1 layered system 
%                  with :1. B1+, B1-, B2+, B2- as boundary condition 
%                  matrices, Phi, Lambda as decomposed layer transfer 
%                  matrix.
%
% [TM] = fun_1layer_pred(B1_pos,B1_neg,B2_pos, B2_neg, Phi, Lambda)
%
% Input Parameters:
%   B1_pos,B1_neg: Boundary condition matrices at upper interface
%   B2_pos,B2_neg: Boundary condition matrices at lower interface
%   Phi, Lambda:   Decomposed tranfer matrix of the last layer
%                  [T] = [Phi][Lambda][Phi]^-1
% Output Parameters:  
%   TM:            2x2 tranfer matrix at one frequencies
% Internal Parameters:
%   A1:             Decomposed global matrix
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
% Email: song520@purdue.edu
% 2022 Fall

% Dimension of the global matrix
N = size(B1_pos,2) + size(B1_neg,2);
B1_h = size(B1_pos,1);

A1 = zeros(N);
A1(1:B1_h,     1:2)    = B1_pos;
A1(1:B1_h,     3:end)  = -B1_neg*Phi;
A1(B1_h+1:end, 3:end)  = B2_pos*Phi/Lambda;
A_inv = inv(A1);

% 2x2 tranfer matrix of the whole system
TM = A_inv(1:2,end-size(B2_neg,1)+1:end)*B2_neg;
end

