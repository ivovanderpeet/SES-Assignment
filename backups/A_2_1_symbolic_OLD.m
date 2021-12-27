clear all; close all; clc;
%% Constants
THin = 240;     % [C] Given
TCout = 110;    % [C] Given
Qdot = 5e7;     % [W] Given
mdotC = 150;    % [kg/s] Given
CpC = 4200;     % [J/kg/K] Assumed constant

CC = mdotC*CpC;
TCin = TCout - Qdot/(mdotC*CpC); % [K] Calculate cold-side inlet temperature

%% Symbolic stuff
syms mdotH CpH CH

CR = CH/CC; % [-] Known that CH < Cc
THout = THin - CC./CH.*(TCout-TCin);

% For counterflow
dT1 = THout - TCin;
dT2 = THin - TCout;
dTlm = (dT2 - dT1)./log(dT2./dT1);

UA = Qdot./dTlm;
NTU = UA./CH; % [-] Known that CH < Cc

% Determine unknowns from NTU == 3
CH =    double(vpasolve(NTU==3, CH));
CR =    double(subs(CR))
THout = double(subs(THout))
dT1 =   double(subs(dT1));
dTlm = double(subs(dTlm));
UA =    double(subs(UA));
NTU =   double(subs(NTU));

% Determine CpH from THout and THin
TmH = (THout + THin)/2;
CpH = XSteam("CpL_T", TmH)*1000

% Calcualte mdotH from CH and CpH
mdotH = CH/CpH

% Effectiveness
eps = (1-exp(-NTU.*(1-CR)))./(1-CR.*exp(-NTU.*(1-CR)));

%% Figures
figure(1)
plot([0,1],[TCout, TCin], [0,1],[THin,THout]);
title('Counterflow configuration')
grid on;
legend('Cold', 'Hot');
xlabel('Length [-]')
ylabel('Temperature [K]')


