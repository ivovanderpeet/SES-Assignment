clear all; close all;
%% Vary hot-side outlet temperature

%% Constants
n = 50;
vec = ones(1,n);

T_H_in = vec*(240 + 273);     % [K] Given
T_C_out = vec*(110 + 273);    % [K] Given
mdot_C = 150;           % [kg/s] Given
Qdot = 50*10^6;         % [J/s] Given

cp_H = 4500;            % [J/K] EngineeringToolbox
cp_C = 4200;            % [J/K]

% T_H_out = 150 + 273;    % [K] Afgeschat
T_H_out = linspace(110, 200, n) + 273;

%% Temperature
T_C_in = T_C_out - Qdot/(mdot_C*cp_C);  % [K] Calculate cold-side inlet temperature

% figure()
% title('Counterflow configuration')
% plot([0,1],[T_C_out, T_C_in], [0,1],[T_H_in,T_H_out]);
% grid on;
% legend('Cold', 'Hot');
% xlabel('Length [-]')
% ylabel('Temperature [K]')

% For counterflow:
dT1 = T_H_out - T_C_in;
dT2 = T_H_in - T_C_out;
dT_lm = (dT2 - dT1)./log(dT2./dT1);

dT_H = -T_H_out + T_H_in;
dT_C = T_C_out - T_C_in;

%% Capacity flows
C_H = Qdot./dT_H;
C_C = mdot_C*cp_C;

Cmin = min(C_C,C_H);
Cmax = max(C_C,C_H);
CR = Cmin./Cmax;

mdot_H = C_H/cp_H;

%% Heat transfer
UA = Qdot./dT_lm;
NTU = UA./Cmin;

eps = (1 - exp(-NTU.*(1-CR)))./(1 - CR.*exp(-NTU.*(1-CR)));


%% Plot
figure()
plot(T_H_out, eps)
title('eps')
grid on

figure()
plot(T_H_out, NTU)
title('NTU')
grid on

figure()
plot(T_H_out, CR)
title('CR')
grid on
