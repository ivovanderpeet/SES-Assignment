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

CR = 0.7; % [-] Choose a value for heat capacity ratio CR NIET HOGER DAN 1

mdot_H = (mdot_C*cp_C)/(cp_H*CR); % [kg/s] Calculate hot-side mass flow rate

T_H_out = T_H_in-CR*(T_C_out-T_C_in);

figure()
title('Counterflow configuration')
plot([0,1],[T_C_out, T_C_in], [0,1],[T_H_in,T_H_out]);
grid on;
legend('Cold', 'Hot');
xlabel('Length [-]')
ylabel('Temperature [K]')

% For counterflow:
dT1 = T_H_out - T_C_in;
dT2 = T_H_in - T_C_out;
dT_lm = (dT2 - dT1)/log(dT2/dT1);


%%
