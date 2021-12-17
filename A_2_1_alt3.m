clear all; close all;
%% Choose mdotH

%% Constants
n = 50;
vec = ones(1,n);
vec2 = ones(1,2*n);

% Temperatures
THin = (240 + 273);     % [K] Given
TCout = (110 + 273);    % [K] Given

% Heat and mass flux
Qdot = 50*10^6;         % [J/s] Given
mdotC = 150;           % [kg/s] Given
cpH = 4500;            % [J/K]
cpC = 4200;            % [J/K]

Cc = mdotC*cpC;

%% Temperature
TCin = TCout - Qdot/(mdotC*cpC); % [K] Calculate cold-side inlet temperature

mdotH = 60;
Ch = mdotH*cpH;
Cmin = min(Cc,Ch);
Cmax = max(Cc,Ch);
CR = Cmin./Cmax

THout = THin - Qdot/Ch % [K] Calculate cold-side inlet temperature
% THout = THin + Cc*(TCout-TCin)/Ch 

% For counterflow:
dT1 = THout - TCin;
dT2 = THin - TCout;
dT_lm = (dT2 - dT1)./log(dT2./dT1);

dTh = -THout + THin;
dTc = TCout - TCin;

%% Heat transfer
UA = Qdot./dT_lm;
NTU = UA./Cmin

eps = (1 - exp(-NTU.*(1-CR)))./(1 - CR.*exp(-NTU.*(1-CR)));

figure()
title('Counterflow configuration')
plot([0,1],[TCout, TCin]-273, [0,1],[THin,THout]-273);
grid on;
legend('Cold', 'Hot');
xlabel('Length [-]')
ylabel('Temperature [K]')
