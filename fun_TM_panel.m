function [Taux,TM] = fun_TM_panel(w, c0, hp, ms, Dp, D, theta)
% fun_TM_panel   returns the transfer matrix a stiff panel
% 
% [Taux, TM] = fun_TM_panel(w, c0, hp, ms, Dp, D, theta)
%
% Input Parameters:
%   w:      Frequency vector [rad/s]
%   c0:     Sound speed [m/s]
%   hp:     Panel thickness [m]
%   ms:     Mass per unit area [kg/m^2]
%   Dp:     Longitudinal stiffness per unit width [Pa m]
%   D:      Flexural stiffness per unit width [Pa m^3]
%   theta:  Incidance angle [degree]
% Output Parameters:  [T] = [Taux] [TM] [Taux]^-1
%   [Taux]:   Unit matrix [I]_{4x4}
%   [TM]:     Transfer matrix
%   
%
% Ref: 
% [1] Song, Guochenhao, Zhuang Mo, and J. Stuart Bolton. "A general and 
%     stable approach to modeling and coupling multilayered acoustical 
%     systems with various types of layers." Journal of Sound and 
%     Vibration 567 (2023): 117898.
%
%
% Written by: 
% Guochenhao Song
% Ray W. Herrick Lab, Purdue University
% Email: guochenhaosong@gmail.com
% 2022 Fall

kx = w/c0*sin(theta);
% Construct the Tau matrix entry-wise
Taux = zeros(4,4,length(w));
TM = zeros(4,4,length(w));

C1 = Dp.*kx.^2 - ms.*w.^2;
C2 = D.*kx.^4 - ms.*w.^2;


% Row 1
TM(1,1,:) = 1;
TM(1,2,:) = -1i*kx*hp;
% Row 2
TM(2,2,:) = 1;
% Row 3
TM(3,1,:) = kx*hp./(2*w).*C1;
TM(3,2,:) = 1./(1i*w).*(kx.^2*hp^2/4.*C1 + C2);
TM(3,3,:) = 1;
TM(3,4,:) = 1i*kx*hp;
% Row 4
TM(4,1,:) = 1./(1i*w).*C1;
TM(4,2,:) = -kx*hp./(2*w).*C1;
TM(4,4,:) = 1;

for count = 1:length(w)
    Taux(:,:,count) = eye(4);
end

end