function [Taux,Lambda, kz] = fun_TM_fluid(w,d, rho, K,theta,c0)
% fun_TM_fluid   returns the decomposed transfer matrix of an fluid layer
% 
% [Taux,Lambda] = fun_TM_fluid(w,d, rho, K,theta,c0)
% [Taux,Lambda, kz] = fun_TM_fluid(w,d, rho, K,theta,c0)
%
% Input Parameters:
%   w:      Frequency vector [rad/s]
%   d:      Thickness [m]
%   rho:    Density [kg/m^3]
%   K:      Fluid bulk modulus [Pa]
%   theta:  Incidance angle [degree]
%   c0:     Sound speed [m/s]
% Output Parameters:  [T] = [Taux] [Lambda] [Taus]^-1
%   [Taus]:   Material properties & incident angle
%   [Lambda]: Wave attenuation matrix
%   [kz]:     Wavenumber kz
%
%   Here, [T]_{fluid} = [Taus][Lambda][Taus]^{-1}
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

% Transverse wavenumber
kt = w/c0*sin(theta);
k = w*sqrt(rho/K);
kz = (k.^2 - kt.^2).^0.5;

% Construct the Tau matrix entry-wise
Taux = zeros(2,2,length(w));
Alpha = zeros(2,length(w));
Lambda = zeros(2,2,length(w));

% Column 1
Taux(1,1,:) = 1;
Taux(2,1,:) = kz./w./rho;
% Column 2
Taux(1,2,:) = 1;
Taux(2,2,:) = -kz./w./rho;

% Wave grow/decay constant
Alpha(1,:) = (-1i.*kz);
Alpha(2,:) = (1i.*kz);

for count = 1:length(w)
    Lambda(:,:,count) = diag(exp(Alpha(:,count)*-d));
end