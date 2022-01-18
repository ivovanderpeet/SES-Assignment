clear all; close all; clc
addpath('C:\Users\ivova\OneDrive - TU Eindhoven\03 Education\Masters Courses\4EM70 Sustainable Energy Sources\SES-Assignment\functions')

% Constants
Cond.DH.Qdot = 50e6; % [W]
Cond.R.Qdot = -Cond.DH.Qdot;

% From well
p(1) = 50;
T(1) = 240;
h(1) = XSteam('h_pT', p(1), T(1));
s(1) = XSteam('s_pT', p(1), T(1));

% After flash tank
p(2) = 18; % [bar] DEZE MAG JE KIEZEN MET STUDENTNUMMER
h(2) = h(1);
T(2) = XSteam('T_ph', p(2), h(2));
s(2) = XSteam('s_ph', p(2), h(2));
Flash.x = XSteam('x_ph', p(2), h(2))

% To turbine (only mass fraction steam)
p(3) = p(2);
T(3) = T(2);
h(3) = XSteam('hV_p', p(3));
s(3) = XSteam('sV_p', p(3));

% From turbine to condensor
p(4) = 2; % [Bar] DEZE MAG JE KIEZEN
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
Cond.R.mdot = Cond.DH.Qdot/Cond.q;
Turb.R.mdot = Cond.R.mdot;
Flash.W.mdot = Cond.R.mdot/Flash.x;

Turb.Qdot = Turb.q*Turb.R.mdot;

%% Condensor
Cond.R.Tin = T(4);
Cond.R.Tout = T(5);
Cond.R.Tm = (Cond.R.Tin+Cond.R.Tout)/2;

Cond.DH.Tout = 110;
Cond.DH.V = 150e-3;  % [m3/s]
Cond.DH.p = p(4);
Cond.DH = getToutCp(Cond.DH);


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