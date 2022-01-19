clear all; close all; clc
addpath('C:\Users\ivova\OneDrive - TU Eindhoven\03 Education\Masters Courses\4EM70 Sustainable Energy Sources\SES-Assignment\functions')

% Constants
Cond.DH.Qdot = 50e6; % [W] DH power
Cond.R.Qdot = -Cond.DH.Qdot;

% From well
p(1) = 50;
T(1) = 240;
h(1) = XSteam('h_pT', p(1), T(1));
s(1) = XSteam('s_pT', p(1), T(1));

% After flash tank
% p(2) = 18.1; % [bar] DEZE MAG JE KIEZEN MET STUDENTNUMMER
p(2) = 12.1;
% p(2) = 30.1;
h(2) = h(1);
T(2) = XSteam('T_ph', p(2), h(2));
s(2) = XSteam('s_ph', p(2), h(2));
Flash.x = XSteam('x_ph', p(2), h(2));
Flash.rhoV = XSteam('rhoV_p',p(2));
Flash.rhoL = XSteam('rhoL_p',p(2));

% To turbine (only mass fraction steam)
p(3) = p(2);
T(3) = T(2);
h(3) = XSteam('hV_p', p(3));
s(3) = XSteam('sV_p', p(3));

% From turbine to condensor
% p(4) = 1.9; % [Bar] DEZE MAG JE KIEZEN MET JE STUDENTNUMMER
p(4) = 5.2;
s(4) = s(3);
h(4) = XSteam('h_ps',p(4),s(4));
T(4) = XSteam('T_ps',p(4),s(4));

% From condensor back to well
p(5) = p(4);
T(5) = T(4);
h(5) = XSteam('hL_p',p(5));
s(5) = XSteam('sL_p',p(5));

% Eenheden fixen
p = p*1e5;
h = h*1e3;
s = s*1e3;

Turb.q = h(3) - h(4);
Cond.q = h(4) - h(5); % [J/kg]

% determine massflows
Cond.R.mdot  = Cond.DH.Qdot/Cond.q;
Flash.W.mdot = Cond.R.mdot/Flash.x;
Turb.R.mdot  = Cond.R.mdot;
Flash.R.mdot = Cond.R.mdot;

Turb.Qdot = Turb.q*Turb.R.mdot;

%% Flash tank
Flash.K = 0.1; % [m/s] from reader
Flash.Vg = Flash.K*sqrt((Flash.rhoL-Flash.rhoV)/Flash.rhoV);
Flash.phig = Flash.R.mdot/Flash.rhoV;
Flash.A = Flash.phig/Flash.Vg;
Flash.D = sqrt(4/pi*Flash.A);
Flash.L = Flash.D*3;
Flash.V = Flash.L*Flash.A;
%% Condensor
Cond.R.Tin = T(4);
Cond.R.Tout = T(5);
Cond.R.Tm = (Cond.R.Tin+Cond.R.Tout)/2;

Cond.DH.Tout = 110;
Cond.DH.V = 150e-3;  % [m3/s]
Cond.DH.p = p(4);
Cond.DH = getToutCp(Cond.DH);

% Determine dTlm and UA
Cond.dT1 = Cond.R.Tout - Cond.DH.Tin;
Cond.dT2 = Cond.R.Tin - Cond.DH.Tout;
Cond.dTlm = (Cond.dT2 - Cond.dT1)./log(Cond.dT2./Cond.dT1);
Cond.UA = -Cond.R.Qdot./Cond.dTlm; 

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
Cond.R.L = 10.3;
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

Cond.Cmin = Cond.DH.mdot*Cond.DH.Cp;
Cond.NTU = Cond.UA/Cond.Cmin;
Cond.CR = 0;
Cond.eps = (1-exp(-Cond.NTU.*(1-Cond.CR)))./(1-Cond.CR.*exp(-Cond.NTU.*(1-Cond.CR)));
Cond.eps2 = (Cond.DH.Tout-Cond.DH.Tin)/(Cond.R.Tin-Cond.DH.Tin);

%%
sat = getSatCurve();

figure(1)
plot([s],[T], 'o-'); hold on
plot(sat.s,sat.T)
title('T-s diagram')
xlabel('Entropy s [kJ/kg/K]')
ylabel('Temperature [C]')
text(s,T,{'1','2','3','4','5'})
grid on

figure(2)
plot(h,T); hold on
plot(sat.h,sat.T);
plot([h(5),h(4)],[Cond.DH.Tin,Cond.DH.Tout]);
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
fprintf('Transfer units:\tNTU = %0.2f\t\t[-]\n',Cond.NTU)
fprintf('Effectiveness:\tepsilon = %0.4f\t\t[-]\n',Cond.eps)
fprintf('Effectiveness2:\tepsilon2 = %0.4f\t\t[-]\n',Cond.eps2)
fprintf('Pres. drop:\t\tdp = %0.2e\t\tdp = %0.2e\t[Pa]\n', Cond.R.dp, Cond.DH.dp)
fprintf('Conv. coef.:\th = %0.2e\t\t\th = %0.2e\t[W/m2/K]\n', Cond.R.h, Cond.DH.h)
fprintf('Nusselt number:\tNu = %0.1f\t\tNu = %0.1f\t\t[-]\n', Cond.R.Nu, Cond.DH.Nu)
fprintf('HX Area:\t\tA = %0.2f\t\tA = %0.2f\t\t[m2]\n', Cond.R.Aht, Cond.DH.Aht)
fprintf('Velocity:\t\tv = %0.2f\t\t\tv = %0.2f\t\t[m/s]\n', Cond.R.v, Cond.DH.v)
fprintf('Reynolds:\t\tRe = %0.3e\t\tRe = %0.3e\t[-]\n\n', Cond.R.Re, Cond.DH.Re)