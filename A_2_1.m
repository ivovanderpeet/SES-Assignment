clear all; close all; clc;
%% A.2.1 Design of a heat-only plant

%% Constants
THin = 240;         % [K] Given
TCout = 110;        % [K] Given
mdotC = 150;        % [kg/s] Given
mdotH = 62.39;    % [kg/s] From symbolic

Qdot = 50*10^6;     % [J/s] Given

CpH = 4.30233e+03;   % [J/K] From symbolic
CpC = 4.18838e+03;   % [J/K] From symbolic

CH = mdotH*CpH;
CC = mdotC*CpC;
Cmin = min(CC,CH);
Cmax = max(CC,CH);
CR = Cmin/Cmax;

%% Calculations
TCin = TCout - Qdot/(mdotC*CpC); % [K] Calculate cold-side inlet temperature
THout = THin - CC/CH*(TCout-TCin);

% For counterflow:
dT1 = THout - TCin;
dT2 = THin - TCout;
dTlm = (dT2 - dT1)/log(dT2/dT1);

% Final parameters
UA = Qdot/dTlm;
NTU = UA/Cmin;

% Check CpH back
TmH = (THout + THin)/2;
XSteam("CpL_T", TmH)*1000;

%% Figure
fprintf('Heat capacities:\tCpH = %0.4e,\tCpC = %0.4e\n', CpH, CpC)
fprintf('Mass flows:\t\t\tmdotH = %3.2f,\t\tmdotC = %3.2f\n', mdotH, mdotC)
fprintf('Nondimensional:\t\tNTU = %0.3f,\t\teps = %0.5f\n', NTU, eps)
fprintf('Heat transfer:\t\tUA = %0.3e\n',UA)

figure()
title('Counterflow configuration')
plot([0,1],[TCout, TCin], [0,1],[THin,THout]);
grid on;
legend('Cold', 'Hot');
xlabel('Length [-]')
ylabel('Temperature [K]')
