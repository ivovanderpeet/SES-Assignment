clear all; close all;
%% A.2.1 Design of a heat-only plant

%% Constants
T_H_in = 240 + 273;     % [K] Given
T_C_out = 110 + 273;    % [K] Given
mdot_C = 150;           % [kg/s] Given
Qdot = 50*10^6;         % [J/s] Given

cp_H = 4771.9;          % [J/K] EngineeringToolbox
cp_C = 4200;            % [J/K]

%% Questions
T_C_in = T_C_out - Qdot/(mdot_C*cp_C); % [K] Calculate cold-side inlet temperature

CR = 0.05:0.05:0.95; % [-] Choose a value for heat capacity ratio CR NIET HOGER DAN 1

mdot_H = (mdot_C*cp_C)./(cp_H*CR); % [kg/s] Calculate hot-side mass flow rate

T_H_out = ones(1,19)*T_H_in-CR*(T_C_out-T_C_in);

% For counterflow:
dT1 = T_H_out - T_C_in;
dT2 = T_H_in - T_C_out;
dT_lm = (dT2 - dT1)./log(dT2./dT1);

%%
UA = Qdot./dT_lm;
NTU = (T_C_out-T_C_in)./dT_lm;
eps = (1-exp(-NTU.*(1-CR)))./(1-CR.*exp(-NTU.*(1-CR)));

%% Plot
figure()
plot(CR, eps)
title('eps')
grid on

figure()
plot(CR, NTU)
title('NTU')
grid on

figure()
plot(CR, mdot_H)
title('mdot_H')
grid on