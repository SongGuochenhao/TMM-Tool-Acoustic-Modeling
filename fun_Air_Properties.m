function [rho0,c0,gamma,eta,Pr,Rg,Dc,m,Sm] = fun_Air_Properties(P0,T,RH)
% fun_Air_Properties
%
% Input Parameters: 
%   P0:     ambient pressure in [Pa], 
%   T:      temperature in [C],
%   RH:     relative humidity
% Output Parameters: 
%   rho0:   air density in [kg/m^3]
%   c0:     sound speed in [m/s]
%   gamma:  air specific heat ratio
%   eta:    air dynamic viscosity in [Pa s]
%   Pr:     air Prandlt number
%   Rg:     air specific gas constant in [J kg^-1 K^-1]
%   Dc:     configurational diffusivity in [m^2/s]
%   m:      air molecule mass in [kg]
%   Sm:     air molecule surface area in [m^2]
%
% Written by: Zhuang Mo
% Updated by: Guochenhao Song, 2022 Fall

%%
% Data and psat https://en.wikipedia.org/wiki/Density_of_air
% Specific gas constant of dry air and wate vapor
Rd = 287.058;
Rv = 461.495;
pv = RH*6.102*10^(7.5*T/(T + 237.8));
pd = P0 - pv;
% Average of dry air and water vapor
Rg = P0/(pd/Rd + pv/Rv);
% rho0 = P0/Rg/(T + 273.15);

% Specific heat ratio of humid air (Wong 1984, JASA)
% A = 5.2e-4 + 4e-5*T + 7.5e-7*T^2 + 4.5e-8*T^3;
% gamma = 1.39984 - A*(RH + 0.125);
% c0 = sqrt(gamma*P0/rho0);

% Standard E 1050
c0 = 20.047*sqrt(273.15+T);
rho0 = 1.290*(P0/101.325)*(273.15/(273.15+T))/1000;
% Sutherland's_law for air viscosity
eta = 1.716e-5*((T+273.15)/273.15)^1.5*(273.15+110.4)/(T+273.15+110.4);
Pr = 0.707;
gamma = 1.4;
% Given air properties from Venegas 2016
Sm = 4.3265e-19;
m = 4.8106e-26;
Dc = 1.35e-10;