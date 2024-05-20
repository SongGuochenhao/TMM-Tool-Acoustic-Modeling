function [Taux,Lambda, K] = fun_TM_solid(w,d,rho,E,nu,theta,c0)
% fun_TM_poro   returns the decomposed transfer matrix of an solid-elastic
%               layer
% 
% [Taux,Lambda] = fun_TM_solid(w,d,rho,E,nu,eta,theta,c0)
% [Taux,Lambda, K] = fun_TM_solid(w,d,rho,E,nu,eta,theta,c0)
%
% Input Parameters:
%   w:      Frequency vector [rad/s]
%   d:      Thickness [m]
%   rho:    Density [kg/m^3]
%   E:      Elastic Modulus [Pa]
%   nu:     Possion's Ratio
%   theta:  Incidance angle [degree]
%   c0:     Sound speed [m/s]
% Output Parameters:  [T] = [Taux] [Lambda] [Taus]^-1
%   [Taus]:   Material properties & incident angle
%   [Lambda]: Wave attenuation matrix
%   [K]:      Two wavenumbers
%
%   Here, [T]_{poro} = [Taus][Lambda][Taus]^{-1}
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

% Transverse wavenumber
k0 = w/c0;
kt = k0*sin(theta);

mu = E/(1+nu)/2;           % Shear modulus, equals to Lame constant miu
lambda = E*nu/(1+nu)/(1-2*nu);   % Lame constant, equals to Lame constant lamda

k1 = w*sqrt(rho/(lambda+2*mu));
k3 = w*sqrt(rho/(mu));
K = [k1; k3]; % wavenumbers

k13 = (k1.^2 - kt.^2).^0.5;
k33 = (k3.^2 - kt.^2).^0.5;

D1 = lambda*(k13.^2 +kt.^2)+2*mu*k13.^2;
D2 = 2*mu*kt;

% Construct the Tau matrix entry-wise
Taux = zeros(4,4,length(w));
Alpha = zeros(4,length(w));
Lambda = zeros(4,4,length(w));

% Column 1
Taux(1,1,:) = w.*kt;
Taux(2,1,:) = w.*k13;
Taux(3,1,:) = -D1;
Taux(4,1,:) = -D2.*k13;
% Column 2
Taux(1,2,:) = w.*kt;
Taux(2,2,:) = -w.*k13;
Taux(3,2,:) = -D1;
Taux(4,2,:) = D2.*k13;
% Column 3
Taux(1,3,:) = -w.*k33;
Taux(2,3,:) = w.*kt;
Taux(3,3,:) = -D2.*k33;
Taux(4,3,:) = D1;
% Column 4
Taux(1,4,:) = w.*k33;
Taux(2,4,:) = w.*kt;
Taux(3,4,:) = D2.*k33;
Taux(4,4,:) = D1;

% Wave grow/decay constant
Alpha(1,:) = (-1i.*k13);
Alpha(2,:) = (1i.*k13);
Alpha(3,:) = (-1i.*k33);
Alpha(4,:) = (1i.*k33);

for count = 1:length(w)
    Lambda(:,:,count) = diag(exp(Alpha(:,count)*-d));
end
end