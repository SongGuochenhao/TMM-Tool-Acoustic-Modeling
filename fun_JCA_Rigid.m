function [rhoeq,Keq] = fun_JCA_Rigid(w,sigma,phi,a1,VCL,TCL,rho0,eta,Pr,gamma,P0)
% fun_JCA_Rigid  returns the equivalent fluid properties corresponds to a 
%                rigid porous medium based on the Johnson-Champoux-Allard 
%                (JCA) Model
%
% Input Parameters:
%   w:      Frequency vector [rad/s]
%   sigma:  Flow resistivity
%   phi:    Porosity
%   a1:     Tortuosity
%   VCL:    Viscous characteristic length [m]
%   TCL:    Thermal characteristic length [m]
%   rho0:   Air (fluid) density [kg/m^3]
%   eta:    Air (fluid) dynamic viscosity [Pa s]
%   Pr:     Air (fluid) Prandlt number
%   gamma:  Air (fluid) specific heat ratio
%   P0:     Ambient pressure [Pa]
% Output Parameters:
%   rhoeq:  Effective fluid density of a rigid-frame porous medium
%   Keq:    Effective bulk modulus of a rigid-frame porous medium
%
% Ref:
% [1] Jean Francois Allard and Noureddine Atalla. Propagation of Sound in 
%     Porous Media: Modelling Sound Absorbing Materials 2e. John Wiley & 
%     Sons, 2009.
% [2] David Linton Johnson, Joel Koplik, and Roger Dashen. Theory of 
%     dynamic permeability and tortuosity in fluid saturated porous media. 
%     Journal of Fluid Mechanics, 176:379–402, 1987.
% [3] Yvan Champoux and Jean F. Allard. Dynamic tortuosity and bulk modulus
%     in air-saturated porous media. Journal of Applied Physics, 
%     70(4):1975–1979, 1991.
%
% Written by: 
% Guochenhao Song
% Ray W. Herrick Lab, Purdue University
% Email: guochenhaosong@gmail.com
% 2022 Fall


% Dynamic Density proposed in[2]
M = sigma*phi/1i./w/a1/rho0;
N = (1 + 4i*eta*rho0*w*a1^2/sigma^2/VCL^2/phi^2).^0.5;
rhoeq = rho0*a1/phi*(1 + M.*N);

% Dynamic Bulk Modulus proposed in [3]
Q = 8i*eta/TCL^2/Pr/rho0./w;
S = (1 + 1i*TCL^2*Pr*rho0*w/16/eta).^0.5;
U = (gamma - 1)./(1 - Q.*S);
Keq = gamma*P0/phi./(gamma - U);
end