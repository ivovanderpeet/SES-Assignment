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
DHcond.Qdot = 50e6; % [W]
Rcond.Qdot = -DHcond.Qdot;

%% Follow the cycle
% 1. Saturated liquid   (Condensor out, pump in)
p(1) = 2;
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

% Convert units
p = p*1e5;
h = h*1e3;
s = s*1e3;

% Heat and work
q(1) = h(1)-h(5);
w(2) = h(2)-h(1);
q(3) = h(3)-h(2);
q(4) = h(4)-h(3);
w(5) = h(5)-h(4);

% Checks after calculating cycle
Recon.mdot = Rcond.Qdot/q(1);
Revap.mdot = Recon.mdot;
Rcond.mdot = Recon.mdot;
x5 = XSteam('x_ps', p(5), s(5));



%% Properties of evaporator

% Rankine cycle-side
Revap.Tin = T(3);
Revap.Tout = T(4);
Revap.Tm = (Revap.Tin+Revap.Tout)/2;
Revap.p = p(3);
Revap.q = q(4);
Revap.Qdot = Revap.q*Revap.mdot;

Revap.Cp = XSteam('Cp_pT', Revap.p/1e5, Revap.Tm*0.99999)*1e3;
Revap.rho = XSteam('rho_pT', Revap.p/1e5, Revap.Tm*0.99999);
Revap.mu = XSteam('my_pT', Revap.p/1e5, Revap.Tm*0.99999);
Revap.k = XSteam('tc_pT', Revap.p/1e5, Revap.Tm*0.99999);

% Well-side
Wevap.mdot = 300; % [kg/s] DEZE MAG JE KIEZEN
Wevap.p = 50e5; % SCHATING MOET NOG VERANDERD WORDEN
Wevap.Tin = 240;
[Revap, Wevap] = getToutCp(Revap,Wevap);

%% Properties of economizer
Wecon.mdot = Wevap.mdot;

% Temperatures 
Recon.Tin = T(2);
Recon.Tout = T(3);
Recon.Tm = (Recon.Tin+Recon.Tout)/2;
Recon.p = p(2);
Recon.q = q(3);
Recon.Qdot = Recon.q*Recon.mdot;

Wecon.Tin = Wevap.Tout;
Wecon.p = 50e5; % SCHATTING MOET NOG VERANDERD WORDEN
[Recon,Wecon] = getToutCp(Recon,Wecon);
% Wecon.Tout = Recon.mdot/(Wecon.mdot*Wecon.Cp)*(h(2)-h(3))+Wecon.Tin;



%% Properties of condenser
% Temperatures
Rcond.Tin = T(5);
Rcond.Tout = T(1);
Rcond.Tm = (Rcond.Tin+Rcond.Tout)/2;
% DHcond.Tin = 30.zoveel
DHcond.Tout = 110;

DHcond.V = 150e-3;  % [m3/s]
DHcond.p = 1.53e5;
DHcond = getDHinfo(DHcond);
%% Figure
sat = getSatCurve();

% figure(1)
% plot([s,s(1)],[T,T(1)]); hold on
% plot(sat.s,sat.T)
% title('T-s diagram')
% xlabel('Entropy s [kJ/kg/K]')
% ylabel('Temperature [C]')
% text(s,T,{'1','2','3','4','5'})
% grid on
% 
figure(2)
plot([h,h(1)],[T,T(1)]); hold on
plot(sat.h,sat.T);
plot([h(1),h(5)],[DHcond.Tin,DHcond.Tout]);
plot([h(2),h(3),h(4)],[Wecon.Tout,Wecon.Tin,Wevap.Tin]);
title('Pinch diagram')
xlabel('Enthalpy h [J/kg]')
ylabel('Temperature [C]')
text(h,T,{'1','2','3','4','5'})
grid on
