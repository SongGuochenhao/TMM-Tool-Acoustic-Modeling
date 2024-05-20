function [B_pos,B_neg] = fun_bc(upper, lower, phi,phi2)
% fun_bc   returns the boundary condition matrix at "+" interface and
%          "-" interface, so that [B_pos][V_pos] = [B_neg][V_neg]
% 
% [B_pos,B_neg] = fun_bc(upper, lower)
% [B_pos,B_neg] = fun_bc(upper, lower, phi)
%
% Input Parameters:
%   upper:  'fluid' / 'rigid wall' / 'solid' / 'stiff panel' / 'poro'
%   lower:  'fluid' / 'rigid wall' / 'solid' / 'stiff panel' / 'poro'
%   phi:    porosity, only necessary if one layer is poro-elastic layer
% Output Parameters:  
%   [B_pos]:   Boundary condition matrix at '+' interface
%   [B_neg]:   Boundary condition matrix at '-' interface
%
% Ref: 
% [1] Song, Guochenhao, Zhuang Mo, and J. Stuart Bolton. "A general and 
%     stable approach to modeling and coupling multilayered acoustical 
%     systems with various types of layers." Journal of Sound and 
%     Vibration 567 (2023): 117898.
%
% Written by: 
% Guochenhao Song
% Ray W. Herrick Lab, Purdue University
% Email: guochenhaosong@gmail.com
% 2022 Fall

% Appendix B.1
if strcmp(upper,'poro') && strcmp(lower, 'fluid')
    B_pos = [0, 0,     0,   1, 0, 0;
             0, 0,     0,   0, 0, 1;
             0, 1-phi, phi, 0, 0, 0;
             0, 0,     0,   0, 1, 0];
    B_neg = [-(1-phi), 0;
             -phi,     0;
             0,        1;
             0,        0];
elseif strcmp(upper,'fluid') && strcmp(lower, 'poro')
    B_neg = [0, 0,     0,   1, 0, 0;
             0, 0,     0,   0, 0, 1;
             0, 1-phi, phi, 0, 0, 0;
             0, 0,     0,   0, 1, 0];
    B_pos = [-(1-phi), 0;
             -phi,     0;
             0,        1;
             0,        0];
% Appendix B.2
elseif strcmp(upper,'poro') && strcmp(lower, 'rigid wall')
    B_pos = [0, 0, 0, 1, 0, 1;
             0, 1, 0, 0, 0, 0;
             0, 0, 1, 0, 0, 0;
             1, 0, 0, 0, 0, 0];
    B_neg = [-1,0;
             0, 0;
             0, 0;
             0, 0];
% Appendix B.3
elseif strcmp(upper,'poro') && strcmp(lower, 'solid')
    B_pos = [1, 0, 0, 0, 0, 0;
             0, 1, 0, 0, 0, 0;
             0, 0, 1, 0, 0, 0;
             0, 0, 0, 1, 0, 1;
             0, 0, 0, 0, 1, 0];
    B_neg = [1, 0, 0, 0;
             0, 1, 0, 0;
             0, 1, 0, 0;
             0, 0, 1, 0;
             0, 0, 0, 1];
elseif strcmp(upper,'solid') && strcmp(lower, 'poro')
    B_neg = [1, 0, 0, 0, 0, 0;
             0, 1, 0, 0, 0, 0;
             0, 0, 1, 0, 0, 0;
             0, 0, 0, 1, 0, 1;
             0, 0, 0, 0, 1, 0];
    B_pos = [1, 0, 0, 0;
             0, 1, 0, 0;
             0, 1, 0, 0;
             0, 0, 1, 0;
             0, 0, 0, 1];
% Appendix B.4
elseif strcmp(upper,'poro') && strcmp(lower, 'stiff panel')
    B_pos = [1, 0, 0, 0, 0, 0;
             0, 1, 0, 0, 0, 0;
             0, 0, 1, 0, 0, 0;
             0, 0, 0, 1, 0, 1;
             0, 0, 0, 0, 1, 0];
    B_neg = [1, 0, 0, 0;
             0, 1, 0, 0;
             0, 1, 0, 0;
             0, 0,-1, 0;
             0, 0, 0, 1];
elseif strcmp(upper,'stiff panel') && strcmp(lower, 'poro')
    B_neg = [1, 0, 0, 0, 0, 0;
             0, 1, 0, 0, 0, 0;
             0, 0, 1, 0, 0, 0;
             0, 0, 0, 1, 0, 1;
             0, 0, 0, 0, 1, 0];
    B_pos = [1, 0, 0, 0;
             0, 1, 0, 0;
             0, 1, 0, 0;
             0, 0,-1, 0;
             0, 0, 0, 1];
% Appendix B.5
elseif strcmp(upper,'solid') && strcmp(lower, 'fluid')
    B_pos = [0, 1, 0, 0;
             0, 0, 1, 0;
             0, 0, 0, 1];
    B_neg = [0, 1;
            -1, 0;
             0, 0];
elseif strcmp(upper,'fluid') && strcmp(lower, 'solid')
    B_neg = [0, 1, 0, 0;
             0, 0, 1, 0;
             0, 0, 0, 1];
    B_pos = [0, 1;
            -1, 0;
             0, 0];   
% Appendix B.6
elseif strcmp(upper,'solid') && strcmp(lower, 'rigid wall')
    B_pos = [1, 0, 0, 0;
             0, 1, 0, 0;
             0, 0, 1, 0];
    B_neg = [0, 0;
             0, 0;
            -1, 0];
% Appendix B.7
elseif strcmp(upper,'solid') && strcmp(lower, 'stiff panel')
    B_pos = [1, 0, 0, 0;
             0, 1, 0, 0;
             0, 0, 1, 0;
             0, 0, 0, 1];
    B_neg = [1, 0, 0, 0;
             0, 1, 0, 0;
             0, 0,-1, 0;
             0, 0, 0, 1];
elseif strcmp(upper,'stiff panel') && strcmp(lower, 'solid')
    B_neg = [1, 0, 0, 0;
             0, 1, 0, 0;
             0, 0, 1, 0;
             0, 0, 0, 1];
    B_pos = [1, 0, 0, 0;
             0, 1, 0, 0;
             0, 0,-1, 0;
             0, 0, 0, 1];
% Appendix B.8
elseif strcmp(upper,'stiff panel') && strcmp(lower, 'fluid')
    B_pos = [0, 1, 0, 0;
             0, 0, 1, 0;
             0, 0, 0, 1];
    B_neg = [0, 1;
             1, 0;
             0, 0];
elseif strcmp(upper,'fluid') && strcmp(lower, 'stiff panel')
    B_neg = [0, 1, 0, 0;
             0, 0, 1, 0;
             0, 0, 0, 1];
    B_pos = [0, 1;
             1, 0;
             0, 0]; 
% Same type layer
elseif strcmp(upper,'poro') && strcmp(lower, 'poro')
    B_pos = eye(6);
    B_neg = eye(6); 
    B_pos(3,2) = 1-phi;     B_neg(3,2) = 1-phi2; 
    B_pos(3,3) = phi;       B_neg(3,3) = phi2; 
    B_pos(4,6) = 1;         B_neg(4,6) = 1; 
    B_pos(6,6) = 1/phi;     B_neg(6,6) = 1/phi2; 
elseif strcmp(upper,'solid') && strcmp(lower, 'solid')
    B_pos = eye(4);
    B_neg = eye(4); 
elseif strcmp(upper,'stiff panel') && strcmp(lower, 'stiff panel')
    B_pos = eye(4);
    B_neg = eye(4); 
elseif strcmp(upper,'fluid') && strcmp(lower, 'fluid')
    B_pos = eye(2);
    B_neg = eye(2); 
else
    error(['Incorrect interface.' newline ...
        'Input should be ''fluid'', ''solid'', ''stiff panel'', ''poro'', or ''rigid wall''']);
end

end

