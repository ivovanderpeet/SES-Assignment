clear all; close all; clc
addpath('C:\Users\ivova\OneDrive - TU Eindhoven\03 Education\Masters Courses\4EM70 Sustainable Energy Sources\SES-Assignment\functions')

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
q = zeros(1,5); % [J/kg]
w = zeros(1,5); % [J/kg]
rho = zeros(1,5);

% Constants
Cond.DH.Qdot = 50e6; % [W]
Cond.R.Qdot = -Cond.DH.Qdot;

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
Econ.R.mdot = Cond.R.Qdot/q(1);
Evap.R.mdot = Econ.R.mdot;
Cond.R.mdot = Econ.R.mdot;
x5 = XSteam('x_ps', p(5)/1e5, s(5)/1e3);

%% Properties of evaporator
% Rankine cycle-side
Evap.R.Tin = T(3);
Evap.R.Tout = T(4);
Evap.R.Tm = (Evap.R.Tin+Evap.R.Tout)/2;
Evap.R.p = p(3);
Evap.R.q = q(4);
Evap.R.Qdot = Evap.R.q*Evap.R.mdot;

Evap.R.Cp = XSteam('Cp_pT', Evap.R.p/1e5, Evap.R.Tm*0.99999)*1e3;
Evap.R.rho = XSteam('rho_pT', Evap.R.p/1e5, Evap.R.Tm*0.99999);
Evap.R.mu = XSteam('my_pT', Evap.R.p/1e5, Evap.R.Tm*0.99999);
Evap.R.k = XSteam('tc_pT', Evap.R.p/1e5, Evap.R.Tm*0.99999);

% Well-side
Evap.W.mdot = 350; % [kg/s] DEZE MAG JE KIEZEN
Evap.W.p = 50e5; % SCHATING MOET NOG VERANDERD WORDEN
Evap.W.Qdot = -Evap.R.Qdot;
Evap.W.Tin = 240;
Evap.W = getToutCp(Evap.W);

%% Properties of economizer
% Temperatures 
Econ.R.Tin = T(2);
Econ.R.Tout = T(3);
Econ.R.Tm = (Econ.R.Tin+Econ.R.Tout)/2;
Econ.R.p = p(2);
Econ.R.q = q(3);
Econ.R.Qdot = Econ.R.q*Econ.R.mdot;

Econ.W.mdot = Evap.W.mdot;
Econ.W.Tin = Evap.W.Tout;
Econ.W.p = 50e5; % SCHATTING MOET NOG VERANDERD WORDEN
Econ.W.Qdot = -Econ.R.Qdot;
Econ.W = getToutCp(Econ.W);

%% Properties of condenser
% Temperatures
Cond.R.Tin = T(5);
Cond.R.Tout = T(1);
Cond.R.Tm = (Cond.R.Tin+Cond.R.Tout)/2;

Cond.DH.Tout = 110;
Cond.DH.V = 150e-3;  % [m3/s]
Cond.DH.p = 1.53e5;
Cond.DH = getToutCp(Cond.DH);

%% Condensor heat exchanger
% Determine dTlm and UA
Cond.dT1 = Cond.R.Tout - Cond.DH.Tin;
Cond.dT2 = Cond.R.Tin - Cond.DH.Tout;
Cond.dTlm = (Cond.dT2 - Cond.dT1)./log(Cond.dT2./Cond.dT1);
Cond.UA = -Cond.R.Qdot./Cond.dTlm; %Cond.R.Qdot = Cond.DH.Qdot

% Determine mean density
Cond.R.rho = nan;
Cond.DH.rho = XSteam('rho_pT', Cond.DH.p/1e5, Cond.DH.Tm);

% Determine average dynamic viscosity
Cond.R.mu = nan;
Cond.DH.mu = XSteam('my_pT', Cond.DH.p/1e5, Cond.DH.Tm);

% Determine thermal conductivity
Cond.R.k = nan;
Cond.DH.k = XSteam('tc_pT', Cond.DH.p/1e5, Cond.DH.Tm);

%Determine heat capacity
Cond.R.Cp = nan;
Cond.DH.Cp = XSteam('Cp_pT', Cond.DH.p/1e5, Cond.DH.Tm)*1000;

% Determine Prandtl number
Cond.R.Pr = nan;
Cond.DH.Pr = Cond.DH.Cp*Cond.DH.mu/Cond.DH.k;

% Volume flow
Cond.R.V = nan;
Cond.DH.V = Cond.DH.mdot/Cond.DH.rho;

%HX
Cond.R.L = 10;
Cond.DH.L = Cond.R.L;
Cond.R.D = 20e-3;
Cond.DH.D = Cond.R.D;

Cond.R.Nplatey = 12;
Cond.DH.Nplatey = Cond.R.Nplatey;
Cond.R.Nplatex = 12;
Cond.DH.Nplatex = Cond.R.Nplatex;

Cond.Wi = Cond.R.Nplatex*Cond.R.D*2;
Cond.Hi = Cond.R.Nplatey*Cond.R.D*2;
Cond.L = Cond.R.L;

Cond.R.rough = 0.015e-3;
Cond.DH.rough = Cond.R.rough;

[Cond.R] = counterDucts(Cond.R);
[Cond.DH] = counterDucts(Cond.DH);

% Heat transfer coef.
Cond.t = 1e-3;
Cond.kwall = 50; % STEEL
Cond.U = (Cond.t/Cond.kwall + 1/Cond.DH.h)^-1;
Cond.A = Cond.UA/Cond.U;

%% Evaporator heat exchanger
% Determine dTlm and UA
Evap.dT1 = Evap.R.Tout - Evap.W.Tin;
Evap.dT2 = Evap.R.Tin - Evap.W.Tout;
Evap.dTlm = (Evap.dT2 - Evap.dT1)./log(Evap.dT2./Evap.dT1);
Evap.UA = -Evap.R.Qdot./Evap.dTlm; %Evap.R.Qdot = Evap.W.Qdot

% Determine mean density
Evap.R.rho = nan;
Evap.W.rho = XSteam('rho_pT', Evap.W.p/1e5, Evap.W.Tm);

% Determine average dynamic viscosity
Evap.R.mu = nan;
Evap.W.mu = XSteam('my_pT', Evap.W.p/1e5, Evap.W.Tm);

% Determine thermal Evapuctivity
Evap.R.k = nan;
Evap.W.k = XSteam('tc_pT', Evap.W.p/1e5, Evap.W.Tm);

%Determine heat capacity
Evap.R.Cp = nan;
Evap.W.Cp = XSteam('Cp_pT', Evap.W.p/1e5, Evap.W.Tm)*1000;

% Determine Prandtl number
Evap.R.Pr = nan;
Evap.W.Pr = Evap.W.Cp*Evap.W.mu/Evap.W.k;

% Volume flow
Evap.R.V = nan;
Evap.W.V = Evap.W.mdot/Evap.W.rho;

%HX
Evap.R.L = 4.2;
Evap.W.L = Evap.R.L;
Evap.R.D = 20e-3;
Evap.W.D = Evap.R.D;

Evap.R.Nplatey = 20;
Evap.W.Nplatey = Evap.R.Nplatey;
Evap.R.Nplatex = 20;
Evap.W.Nplatex = Evap.R.Nplatex;

Evap.Wi = Evap.R.Nplatex*Evap.R.D*2;
Evap.Hi = Evap.R.Nplatey*Evap.R.D*2;
Evap.L = Evap.R.L;

Evap.R.rough = 0.015e-3;
Evap.W.rough = Evap.R.rough;

[Evap.R] = counterDucts(Evap.R);
[Evap.W] = counterDucts(Evap.W);

% Heat transfer coef.
Evap.t = 1e-3;
Evap.kwall = 50; % STEEL
Evap.U = (Evap.t/Evap.kwall + 1/Evap.W.h)^-1;
Evap.A = Evap.UA/Evap.U;

%% Economizer heat exchanger
% Determine dTlm and UA
Econ.dT1 = Econ.R.Tout - Econ.W.Tin;
Econ.dT2 = Econ.R.Tin - Econ.W.Tout;
Econ.dTlm = (Econ.dT2 - Econ.dT1)./log(Econ.dT2./Econ.dT1);
Econ.UA = -Econ.R.Qdot./Econ.dTlm; %Econ.R.Qdot = Econ.W.Qdot

% Determine mean density
Econ.R.rho = XSteam('rho_pT', Econ.R.p/1e5, Econ.R.Tm);
Econ.W.rho = XSteam('rho_pT', Econ.W.p/1e5, Econ.W.Tm);

% Determine average dynamic viscosity
Econ.R.mu = XSteam('my_pT', Econ.R.p/1e5, Econ.R.Tm);
Econ.W.mu = XSteam('my_pT', Econ.W.p/1e5, Econ.W.Tm);

% Determine thermal Econuctivity
Econ.R.k = XSteam('tc_pT', Econ.R.p/1e5, Econ.R.Tm);
Econ.W.k = XSteam('tc_pT', Econ.W.p/1e5, Econ.W.Tm);

%Determine heat capacity
Econ.R.Cp = XSteam('Cp_pT', Econ.R.p/1e5, Econ.R.Tm)*1000;
Econ.W.Cp = XSteam('Cp_pT', Econ.W.p/1e5, Econ.W.Tm)*1000;

% Determine Prandtl number
Econ.R.Pr = Econ.R.Cp*Econ.R.mu/Econ.R.k;
Econ.W.Pr = Econ.W.Cp*Econ.W.mu/Econ.W.k;

% Volume flow
Econ.R.V = Econ.R.mdot/Econ.R.rho;
Econ.W.V = Econ.W.mdot/Econ.W.rho;

% Counterflow ducts HX
% Econ.R.L = 3.7;
% Econ.W.L = Econ.R.L;
% Econ.R.D = 15e-3;
% Econ.W.D = Econ.R.D;
% 
% Econ.R.Nplatey = 10;
% Econ.W.Nplatey = Econ.R.Nplatey;
% Econ.R.Nplatex = 10;
% Econ.W.Nplatex = Econ.R.Nplatex;
% 
% Econ.Wi = Econ.R.Nplatex*Econ.R.D*2;
% Econ.Hi = Econ.R.Nplatey*Econ.R.D*2;
% Econ.L = Econ.R.L;
% 
% Econ.R.rough = 0.015e-3;
% Econ.W.rough = Econ.R.rough;
% 
% [Econ.R] = counterDucts(Econ.R);
% [Econ.W] = counterDucts(Econ.W);

% Heat transfer coef.
% Econ.t = 1e-3;
% Econ.kwall = 50; % STEEL
% Econ.U = (1/Econ.R.h + Econ.t/Econ.kwall + 1/Econ.W.h)^-1;
% Econ.A = Econ.UA/Econ.U;

% Shell & tube HX
Econ.R.ID = 20e-3;
Econ.R.OD = Econ.R.ID + 2* 1e-3;
Econ.R.pitch = 1.5*Econ.R.OD; % Distance between centers of two tubes 1.25 OD is conventional

Econ.R.Ntube = 75;
Econ.R.Npass = 6;
Econ.R.L = 2;

Econ.R.rough = 0.015e-3;
Econ.W.rough = Econ.R.rough;

% https://www.engineeringtoolbox.com/smaller-circles-in-larger-circle-d_1849.html
% http://hydra.nat.uni-magdeburg.de/packing/cci/
load('circ.mat');
Econ.W.ID = Econ.R.pitch*circ.ratio(Econ.R.Ntube*Econ.R.Npass);
Econ.W.L = Econ.R.L;

[Econ.R,Econ.W] = shellTube(Econ.R,Econ.W);

Econ.kWall = 50; % Steel = 50 W/m/K
Econ.U = (1/Econ.W.h + Econ.R.OD*log(Econ.R.OD/Econ.R.ID)/2/Econ.kWall + Econ.R.OD/Econ.R.ID/Econ.R.h)^-1;  % 4PC00 7.43
Econ.A = Econ.UA/Econ.U;

%% Kostenplaatje
costTurbCond = 1e6;
% costCompress = 
% cost = 



%% Figure
sat = getSatCurve();
clc;

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
plot([h(1),h(5)],[Cond.DH.Tin,Cond.DH.Tout]);
plot([h(2),h(3),h(4)],[Econ.W.Tout,Econ.W.Tin,Evap.W.Tin]);
title('Pinch diagram')
xlabel('Enthalpy h [J/kg]')
ylabel('Temperature [C]')
text(h,T,{'1','2','3','4','5'})
grid on

fprintf('CONDENSOR HEAT EXCHANGER PROPERTIES\n')
fprintf('Heat-trns coef:\tU = %0.2f\t\t[W/m2/K]\n',Cond.U)
fprintf('Reqired area:\tA = %0.2f\t\t[m2]\n',Cond.A)
fprintf('Total length:\tL = %0.2f\t\t[m]\n',Cond.L)
fprintf('Total width:\tW = %0.2f\t\t[m]\n',Cond.Wi)
fprintf('Total height:\tH = %0.2f\t\t[m]\n',Cond.Hi)
fprintf('Pres. drop:\t\tdp = %0.2e\t\tdp = %0.2e\t[Pa]\n', Cond.R.dp, Cond.DH.dp)
fprintf('Conv. coef.:\th = %0.2e\t\t\th = %0.2e\t[W/m2/K]\n', Cond.R.h, Cond.DH.h)
fprintf('Nusselt number:\tNu = %0.1f\t\tNu = %0.1f\t\t[-]\n', Cond.R.Nu, Cond.DH.Nu)
fprintf('HX Area:\t\tA = %0.2f\t\tA = %0.2f\t\t[m2]\n', Cond.R.Aht, Cond.DH.Aht)
fprintf('Velocity:\t\tv = %0.2f\t\t\tv = %0.2f\t\t[m/s]\n', Cond.R.v, Cond.DH.v)
fprintf('Reynolds:\t\tRe = %0.3e\t\tRe = %0.3e\t[-]\n\n', Cond.R.Re, Cond.DH.Re)

fprintf('EVAPORATOR HEAT EXCHANGER PROPERTIES\n')
fprintf('Heat-trns coef:\tU = %0.2f\t[W/m2/K]\n',Evap.U)
fprintf('Reqired area:\tA = %0.2f\t\t[m2]\n',Evap.A)
fprintf('Total length:\tL = %0.2f\t\t[m]\n',Evap.L)
fprintf('Total width:\tW = %0.2f\t\t[m]\n',Evap.Wi)
fprintf('Total height:\tH = %0.2f\t\t[m]\n',Evap.Hi)
fprintf('Pres. drop:\t\tdp = %0.2e\tdp = %0.2e\t[Pa]\n', Evap.W.dp, Evap.R.dp)
fprintf('Conv. coef.:\th = %0.2e\th = %0.2e\t[W/m2/K]\n', Evap.W.h, Evap.R.h)
fprintf('Nusselt number:\tNu = %0.1f\t\tNu = %0.1f\t\t[-]\n', Evap.W.Nu, Evap.R.Nu)
fprintf('HX Area:\t\tA = %0.2f\t\tA = %0.2f\t\t[m2]\n', Evap.W.Aht, Evap.R.Aht)
fprintf('Velocity:\t\tv = %0.2f\t\tv = %0.2f\t\t[m/s]\n', Evap.W.v, Evap.R.v)
fprintf('Reynolds:\t\tRe = %0.3e\tRe = %0.3e\t[-]\n\n', Evap.W.Re, Evap.R.Re)

fprintf('ECONOMIZER HEAT EXCHANGER PROPERTIES\n')
fprintf('Heat-trns coef:\tU = %0.2f\t\t[W/m2/K]\n',Econ.U)
fprintf('Reqired area:\tA = %0.2f\t\t[m2]\n',Econ.A)
fprintf('Total length:\tL = %0.2f\t\t[m]\n',Econ.R.L)
fprintf('Shell diameter:\tD = %0.2f\t\t[m]\n',Econ.W.ID)
fprintf('Tube diameter:\tID = %0.4f\t\t[m]\n',Econ.R.ID)
fprintf('Tube diameter:\tOD = %0.4f\t\t[m]\n',Econ.R.OD)
fprintf('Num. of tubes:\tN = %0.2f\t\t[-]\n',Econ.R.Ntube)
fprintf('Num. of passes:\tN = %0.2f\t\t[-]\n',Econ.R.Npass)
fprintf('Pres. drop:\t\tdp = %0.2e\tdp = %0.2e\t[Pa]\n', Econ.W.dp, Econ.R.dp)
fprintf('Conv. coef.:\th = %0.2e\th = %0.2e\t[W/m2/K]\n', Econ.W.h, Econ.R.h)
fprintf('Nusselt number:\tNu = %0.1f\t\tNu = %0.1f\t\t[-]\n', Econ.W.Nu, Econ.R.Nu)
fprintf('HX Area:\t\tA = %0.2f\t\tA = %0.2f\t\t[m2]\n', Econ.W.Aht, Econ.R.Aht)
fprintf('Velocity:\t\tv = %0.2f\t\tv = %0.2f\t\t[m/s]\n', Econ.W.v, Econ.R.v)
fprintf('Reynolds:\t\tRe = %0.3e\tRe = %0.3e\t[-]\n\n', Econ.W.Re, Econ.R.Re)
