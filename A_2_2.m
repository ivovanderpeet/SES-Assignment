clear all; close all; clc
% 1. Saturated liquid   (Condensor out, pump in)
% 2. Compressed liquid  (Pump out, economizer in)
% 3. Saturated liquid   (Economizer out, evaporator in)
% 4. Saturated vapour   (Evaporator out, turbine in)
% 5. Liquid & vapour    (Turbine out, condensor in)

%% Initializing
p = zeros(1,5); % [bar]
T = zeros(1,5); % [C]
h = zeros(1,5); % [kJ/kg]
s = zeros(1,5); % [kJ/kg/K]
q = zeros(1,5); % [kJ/kg]
w = zeros(1,5); % [kJ/kg]
rho = zeros(1,5);

% Constants
Qdot = 50e3; % [kW]

%% Follow the cycle
% 1. Saturated liquid   (Condensor out, pump in)
p(1) = 1.5;
T(1) = XSteam('Tsat_p', p(1));
h(1) = XSteam('hL_p', p(1));
s(1) = XSteam('sL_p', p(1));
rho(1) = XSteam('rhoL_p', p(1));

% 3. Saturated liquid   (Economizer out, evaporator in)
T(3) = 200;
p(3) = XSteam('psat_T', T(3));
h(3) = XSteam('hL_T', T(3));
s(3) = XSteam('sL_T', T(3));
rho(3) = XSteam('rhoL_p', p(3));

% 4. Saturated vapour   (Evaporator out, turbine in)
T(4) = T(3);
p(4) = XSteam('psat_T', T(4));
h(4) = XSteam('hV_T', T(4));
s(4) = XSteam('sV_T', T(4));
rho(4) = XSteam('rhoV_p', p(4));

% 5. Liquid & vapour    (Turbine out, condensor in)
p(5) = p(1);
T(5) = T(1);
s(5) = s(4);
h(5) = XSteam('h_ps', p(5), s(5));
rho(5) = XSteam('rho_ps', p(5), s(5));

% 2. Compressed liquid  (Pump out, economizer in)
p(2) = p(3);
s(2) = s(1);
T(2) = XSteam('T_ps', p(2), s(2));
h(2) = XSteam('h_ps', p(2), s(2));
rho(2) = XSteam('rho_ps', p(2), s(2));

% Heat and work
q(1) = h(1)-h(5);
w(2) = h(2)-h(1);
q(3) = h(3)-h(2);
q(4) = h(4)-h(3);
w(5) = h(5)-h(4);

%% Checks after calculating cycle
mdot = -Qdot/q(1);
x5 = XSteam('x_ps', p(5), s(5));

%% Figure
sat = getSatCurve();

figure(1)
plot([s,s(1)],[T,T(1)]); hold on
plot(sat.s,sat.T)
title('T-s diagram')
xlabel('Entropy s [kJ/kg/K]')
ylabel('Temperature [C]')
text(s,T,{'1','2','3','4','5'})
grid on

figure(2)
plot([h,h(1)],[T,T(1)]); hold on
plot(sat.h,sat.T)
title('Pinch diagram')
xlabel('Enthalpy h [J/kg]')
ylabel('Temperature [C]')
text(h,T,{'1','2','3','4','5'})
grid on
