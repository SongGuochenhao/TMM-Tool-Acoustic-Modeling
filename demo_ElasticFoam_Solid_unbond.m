% demo_ElasticFoam_Solid_unbond    Predict the absorption and transmission
%                                  spectrum of a double layere system:
%
%          Air
%    ------------------
%    |                |
%    |  30 mm foam    |  Top layer
%    |                |
%    ------------------  Unbonded interface
%    |  4 mm plastic  |  Bottom layer
%    ------------------
%          Air
%
% Note that this is an example case in the following reference:
% Song, Guochenhao, Zhuang Mo, and J. Stuart Bolton. "A general and stable 
% approach to modeling and coupling multilayered acoustical systems with 
% various types of layers." Journal of Sound and Vibration 567 (2023): 117898.
%
% Written by: 
% Guochenhao Song
% Ray W. Herrick Lab, Purdue University
% Email: guochenhaosong@gmail.com
% 2022 Fall

clc; clear all

% Frequency in interest
f_idx = linspace(log10(100),log10(6.4e3),1000);
f = 10.^f_idx;
w = 2*pi*f;

% Incident angle
theta = 50/180*pi;


% Air properties
P0 = 101325;
T = 20;
RH = 0.2;
[rho0,c0,gamma,eta,Pr] = fun_Air_Properties(P0, T, RH);
z0 = rho0*c0;


% material properties
% melamine foam
foam_h = 30e-3;
foam_VCL = 199e-6;
foam_TCL = 445e-6;
foam_TOR = 1.00;
foam_phi = 0.99;
foam_rho1 = 10.3;
foam_sigma = 1.31e4;
foam_eta = 0.08;
foam_E = 1.32e5*(1+1j*foam_eta);
foam_nu = 0.33;

% plastic
PP_h = 4e-3;
PP_rho = 920;
PP_E = 1.3e9;
PP_est = 0;
PP_nu = 0.43;

% BC matrices
[B1_pos, B1_neg] = fun_bc('fluid','poro',foam_phi);
[B2_pos, B2_neg] = fun_bc('poro','fluid',foam_phi);
[B3_pos, B3_neg] = fun_bc('fluid','solid');
[B4_pos, B4_neg] = fun_bc('solid','fluid');

% for each freq
TM = zeros(2,2,length(w));
tc = zeros(1,length(w));
rc = zeros(1,length(w));
for count = 1:length(w)
    % upper layer - poro-elastic foam 
    [rhoeq,Keq] = fun_JCA_Rigid(w(count),foam_sigma,foam_phi,foam_TOR,...
    foam_VCL,foam_TCL,rho0,eta,Pr,gamma,P0);
    sigma = foam_sigma;
    rho22 = foam_phi^2*rhoeq;
    rho12 = foam_phi*rho0 - rho22;
    rho11 = foam_rho1 - rho12;
    [Phi,Lambda] =  fun_TM_poro(w(count),foam_h,...
            foam_phi,foam_E,Keq,foam_nu,theta,c0,rho11,rho12,rho22);
    % simplify upper layer to 2x2 TM
    TM1 = fun_1layer_pred(B1_pos,B1_neg,B2_pos, B2_neg, Phi, Lambda);
    % lower layer - elastic-solid
    [Phi,Lambda] = fun_TM_solid(w(count),PP_h,PP_rho,...
        PP_E,PP_nu,theta,c0);
    % simplify lower layer to 2x2 TM
    TM2 = fun_1layer_pred(B3_pos,B3_neg,B4_pos, B4_neg, Phi, Lambda);

    % 2x2 TM for the whole system
    TM(:,:,count) = TM1*TM2;
    d = PP_h+foam_h;
    
    % Solve for transmission & reflection based on 2x2 transfer matrix
    Denom = TM(1,1,count)+TM(1,2,count)*cos(theta)/z0+TM(2,1,count)*z0/cos(theta)+TM(2,2,count);
    tc(count) = 2*exp(1i*w(count)*cos(theta)/c0*d)./Denom;
    rc(count) = (TM(1,1,count)+TM(1,2,count)*cos(theta)/z0-TM(2,1,count)*z0/cos(theta)-TM(2,2,count))./Denom;
end
% absorption and transmission spectra
alpha = 1 - abs(rc).^2;
TL = 20*log10(1./abs(tc));

% plot absorption & transmission
subplot 211
semilogx(f,alpha,'linewidth',2);
xlim(round([f(1),f(end)]))
xticks([100,200,400,800,1600,3200,6400])
hold on
xlabel('Frequency - Hz');
ylabel('Absorption coefficient')
ylim([0 1]);
grid on;
set(0,'DefaultAxesFontName', 'Times New Roman');
set(gca,'FontSize',22);
set(gcf,'color','white')
yticks([0:.2:1])
subplot 212
semilogx(f,TL,'linewidth',2);
hold on
xlabel('Frequency - Hz');
ylabel('Transmission loss - dB');
grid on;
yticks([0:20:200]);
xticks([100,200,400,800,1600,3200,6400])
xlim(round([f(1),f(end)]))
set(0,'DefaultAxesFontName', 'Times New Roman');
set(gca,'FontSize',22);