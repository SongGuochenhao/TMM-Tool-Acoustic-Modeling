function [Taux,Lambda,K] = fun_TM_poro(w,d,phi, E, Kf,nu,theta,c0,rho11,rho12,rho22)
% fun_TM_poro   returns the decomposed transfer matrix of an poro-elastic
%               layer
% 
% [Taux,Lambda] = fun_TM_poro(w,d,phi, E, Kf,nu,theta,c0,rho11,rho12,rho22)
% [Taux,Lambda,K] = fun_TM_poro(w,d,phi, E, Kf,nu,theta,c0,rho11,rho12,rho22)
%
% Input Parameters:
%   w:      Frequency vector [rad/s]
%   d:      Thickness [m]
%   phi:    Open porosity
%   E:      Elastic Modulus [Pa]
%   Kf:     Fluid bulk modulus [Pa]
%   nu:     Possion's Ratio
%   theta:  Incidance angle [degree]
%   c0:     Sound speed [m/s]
%   rho11, rho12, rho22:  Density [kg/m^3]
% Output Parameters:  [T] = [Taux] [Lambda] [Taus]^-1
%   [Taus]:   Material properties & incident angle
%   [Lambda]: Wave attenuation matrix
%   [K]:      Three wavenumbers
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
kt = w/c0*sin(theta);
Kb = E/3/(1-2*nu);
N = E/2/(1 + nu);
% Kb = 2*N*(1 + nu)/3/(1 - 2*nu);
P = 4/3*N + Kb + Kf*(1-phi).^2;
Q = Kf*phi*(1-phi);
R = Kf*phi^2;
Delta = (P.*rho22 + R.*rho11-2*Q.*rho12).^2 - 4*(P.*R-Q.^2).*(rho11.*rho22 - rho12.^2);
k1 = ((w.^2.*(P.*rho22 + R.*rho11-2*Q.*rho12 - Delta.^0.5))/2./(P.*R-Q.^2)).^0.5;
k2 = ((w.^2.*(P.*rho22 + R.*rho11-2*Q.*rho12 + Delta.^0.5))/2./(P.*R-Q.^2)).^0.5;
k3 = (w.^2.*(rho11.*rho22 - rho12.^2)/N./rho22).^0.5;
K = [k1; k2; k3]; % wavenumbers
k13 = (k1.^2 - kt.^2).^0.5;
k23 = (k2.^2 - kt.^2).^0.5;
k33 = (k3.^2 - kt.^2).^0.5;
mu1 = (P.*k1.^2-w.^2.*rho11)./(w.^2.*rho12-Q.*k1.^2);
mu2 = (P.*k2.^2-w.^2.*rho11)./(w.^2.*rho12-Q.*k2.^2);
mu3 = -rho12./rho22;
D1 = (P + Q.*mu1).*k1.^2 - 2*N.*kt.^2;
D2 = (P + Q.*mu2).*k2.^2 - 2*N.*kt.^2;
E1 = (R.*mu1 + Q).*k1.^2;
E2 = (R.*mu2 + Q).*k2.^2;

% Construct the Tau matrix entry-wise
Taux = zeros(6,6,length(w));
Alpha = zeros(6,length(w));
Lambda = zeros(6,6,length(w));

% Column 1
Taux(1,1,:) = w.*kt;
Taux(2,1,:) = w.*k13;
Taux(3,1,:) = w.*k13.*mu1;
Taux(4,1,:) = -D1;
Taux(5,1,:) = -2*N.*kt.*k13;
Taux(6,1,:) = -E1;
% Column 2
Taux(1,2,:) = w.*kt;
Taux(2,2,:) = -w.*k13;
Taux(3,2,:) = -w.*mu1.*k13;
Taux(4,2,:) = -D1;
Taux(5,2,:) = 2*N.*kt.*k13;
Taux(6,2,:) = -E1;
% Column 3
Taux(1,3,:) = w.*kt;
Taux(2,3,:) = w.*k23;
Taux(3,3,:) = w.*k23.*mu2;
Taux(4,3,:) = -D2;
Taux(5,3,:) = -2*N.*kt.*k23;
Taux(6,3,:) = -E2;
% Column 4
Taux(1,4,:) = w.*kt;
Taux(2,4,:) = -w.*k23;
Taux(3,4,:) = -w.*mu2.*k23;
Taux(4,4,:) = -D2;
Taux(5,4,:) = 2*N.*kt.*k23;
Taux(6,4,:) = -E2;
% Column 5
Taux(1,5,:) = -w.*k33;
Taux(2,5,:) = w.*kt;
Taux(3,5,:) = w.*kt.*mu3;
Taux(4,5,:) = -2*N.*k33.*kt;
Taux(5,5,:) = N.*(k33.^2 - kt.^2);
% Column 6
Taux(1,6,:) = w.*k33;
Taux(2,6,:) = w.*kt;
Taux(3,6,:) = w.*kt.*mu3;
Taux(4,6,:) = 2*N.*k33.*kt;
Taux(5,6,:) = N.*(k33.^2 - kt.^2);

% Wave grow/decay constant
Alpha(1,:) = (-1i.*k13);
Alpha(2,:) = (1i.*k13);
Alpha(3,:) = (-1i.*k23);
Alpha(4,:) = (1i.*k23);
Alpha(5,:) = (-1i.*k33);
Alpha(6,:) = (1i.*k33);

for count = 1:length(w)
    Lambda(:,:,count) = diag(exp(Alpha(:,count)*-d));
end