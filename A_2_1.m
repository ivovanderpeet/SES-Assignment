clear all; close all; clc;
addpath('C:\Users\ivova\OneDrive - TU Eindhoven\03 Education\Masters Courses\4EM70 Sustainable Energy Sources\SES-Assignment\functions')


%% Constants
H.Tin = 240;     % [C] Given
C.Tout = 110;    % [C] Given
C.Qdot = 5e7;      % [W] Given
C.V = 150e-3;    % [kg/s] Given
H.p = 50e5;
C.p = 1.53e5;

%% Loop for determining C.Cp
C = getToutCp(C);

%% Symbolic hot-side stuff
syms CH
H.C = CH;

CR = C.C/H.C; % [-] Given that C.C < H.C
H.Tout = H.Tin - C.C./H.C.*(C.Tout-C.Tin);

% For counterflow
dT1 = H.Tout - C.Tin;
dT2 = H.Tin - C.Tout;
dTlm = (dT2 - dT1)./log(dT2./dT1);

UA = C.Qdot./dTlm;
NTU = UA./C.C; % [-] Given that C.C < H.C

% Determine unknowns from CR = 1
H.C = double(vpasolve(CR==1, H.C));
CH = H.C;
CR =    double(subs(CR));
H.Tout = double(subs(H.Tout));
dT1 =   double(subs(dT1));

if CR == 1
    dTlm = dT1;
    UA = C.Qdot/dTlm;
    NTU = UA/C.C;
else
    dTlm =  double(subs(dTlm));
    UA =    double(subs(UA));
    NTU =   double(subs(NTU));
end

clear('CH');

% Determine H.Cp from H.Tout and H.Tin
H.Tm = (H.Tout + H.Tin)/2;
H.Cp = XSteam('Cp_pT', H.p/1e5, H.Tm)*1000;

% Calcualte H.mdot from H.C and H.Cp
H.mdot = H.C/H.Cp;

% Calculate effectiveness
if CR == 1
    eps = NTU/(NTU+1);
else
    eps = (1-exp(-NTU.*(1-CR)))./(1-CR.*exp(-NTU.*(1-CR)));
end

% Determine mean density
H.rho = XSteam('rho_pT', H.p/1e5, H.Tm);
C.rho = XSteam('rho_pT', C.p/1e5, C.Tm);

% Determine average dynamic viscosity
H.mu = XSteam('my_pT', H.p/1e5, H.Tm);
C.mu = XSteam('my_pT', C.p/1e5, C.Tm);

% Determine thermal conductivity
H.k = XSteam('tc_pT', H.p/1e5, H.Tm);
C.k = XSteam('tc_pT', C.p/1e5, C.Tm);

% Determine Prandtl number
H.Pr = H.Cp*H.mu/H.k;
C.Pr = C.Cp*C.mu/C.k;

% Volume flow
H.V = H.mdot/H.rho;
C.V = C.mdot/C.rho;

%% Shell & tube heat exchanger
H.ST.ID = 20e-3;
H.ST.OD = H.ST.ID + 2* 2e-3;
H.ST.pitch = 1.25*H.ST.OD; % Distance between centers of two tubes 1.25 OD is conventional
H.ST.Ntube = 210;
H.ST.Npass = 2;
H.ST.L = 5;
H.ST.rough = 0.015e-3;
C.ST.rough = H.ST.rough;

% https://www.engineeringtoolbox.com/smaller-circles-in-larger-circle-d_1849.html
% http://hydra.nat.uni-magdeburg.de/packing/cci/
load('circ.mat');
C.ST.ID = H.ST.pitch*circ.ratio(H.ST.Ntube*H.ST.Npass);
C.ST.L = H.ST.L;

% Unsure about C.h (effective diameter)
[H,C] = shellTube(H,C);

% Correction factor Shell & Tube
ST.Z = (H.Tin-H.Tout)/(C.Tout-C.Tin);
ST.Y = (H.Tout-H.Tin)/(C.Tin-H.Tin);
ST.F = 0.93;

ST.kWall = 50; % Steel = 50 W/m/K
ST.U = (1/C.ST.h + H.ST.OD*log(H.ST.OD/H.ST.ID)/2/ST.kWall + H.ST.OD/H.ST.ID/H.ST.h)^-1;  % 4PC00 7.43
ST.A = UA/ST.U;
ST.AF = ST.A/ST.F;

%% Counterflow plate heat exchanger
H.CF.L = 3;
C.CF.L = H.CF.L;
H.CF.H = 10e-3;
C.CF.H = H.CF.H;
H.CF.W = 0.5;
C.CF.W = H.CF.W;
H.CF.Nplate = 20;
C.CF.Nplate = H.CF.Nplate;

CF.H = H.CF.H*H.CF.Nplate*2;
CF.W = H.CF.W;
CF.L = H.CF.L;

H.CF.rough = 0.015e-3;
C.CF.rough = H.CF.rough;

[H,C] = counterPlate(H,C);

% Heat transfer coef.
CF.t = 2e-3;
CF.kwall = 50; % STEEL
CF.U = (1/H.CF.h + CF.t/CF.kwall + 1/C.CF.h)^-1;
CF.A = UA/CF.U;

%% Counterflow ducts heat exchanger
CD.C = C;
CD.H = H;
CD.H.L = 3.7;
CD.C.L = CD.H.L;
CD.H.D = 20e-3;
CD.C.D = CD.H.D;

CD.H.Nplatey = 15;
CD.C.Nplatey = CD.H.Nplatey;
CD.H.Nplatex = 15;
CD.C.Nplatex = CD.H.Nplatex;

CD.W = CD.H.Nplatex*CD.H.D*2;
CD.Height = CD.H.Nplatey*CD.H.D*2;
CD.L = CD.H.L;

CD.H.rough = 0.015e-3;
CD.C.rough = CD.H.rough;

[CD.H,CD.C] = counterDucts(CD.H,CD.C);

% Heat transfer coef.
CD.t = 1e-3;
CD.kwall = 50; % STEEL
CD.U = (1/CD.H.h + CD.t/CD.kwall + 1/CD.C.h)^-1;
CD.A = UA/CD.U;

%% Validate calculations
eps2 = (C.Tout-C.Tin)/(H.Tin-C.Tin);


%% Figures & display
clc
% fprintf('THERMAL PROPERTIES OF FLUIDS\n')
fprintf('\t\t\t\tHOT SIDE\t\tCOLD SIDE\t\tUNIT\n')
fprintf('Temp. in:\t\tTin = %0.2f\tTin = %0.2f\t\t[deg C]\n', H.Tin, C.Tin)
fprintf('Temp. out:\t\tTout = %0.2f\tTout = %0.2f\t[deg C]\n', H.Tout, C.Tout)
fprintf('Mean temp.:\t\tTm = %0.2f\t\tTm = %0.2f\t\t[deg C]\n\n', H.Tm, C.Tm)
fprintf('\t\t\t\tHOT SIDE\t\tCOLD SIDE\t\tUNIT\n')
fprintf('Heat capacity:\tCp = %0.2f\tCp = %0.2f\t[J/kg/K]\n', H.Cp, C.Cp)
fprintf('Mass flows:\t\tmdot = %3.2f\tmdot = %3.2f\t[kg/s]\n', H.mdot, C.mdot)
fprintf('Liq. density:\trho = %0.2f\trho = %0.2f\t[kg/m3]\n',H.rho, C.rho)
fprintf('Dyn. visc.:\t\tmu = %0.3e\tmu = %0.3e\t[Pa s]\n',H.mu, C.mu)
fprintf('Therm. cond.:\tk = %0.4f\t\tk = %0.4f\t\t[W/m/K]\n\n',H.k, C.k)

fprintf('OVERALL & NONDIMENSIONAL PROPERTIES\n')
fprintf('Heat transfer:\tUA = %0.3e\t[W/K]\n',UA)
fprintf('Transfer units:\tNTU = %0.2f\t\t[-]\n',NTU)
fprintf('Effectiveness:\teps = %0.4f\t[-]\n', eps)
fprintf('Effectiveness2:\teps2 = %0.4f\t[-]\n\n', eps2)

fprintf('SHELL & TUBE HEAT EXCHANGER PROPERTIES\n')
fprintf('Corr. factor:\tZ = %0.3f\t\t[-]\n',ST.Z)
fprintf('Corr. factor:\tY = %0.3f\t\t[-]\n',ST.Y)
fprintf('Corr. factor:\tF = %0.3f\t\t[-]\n',ST.F)
fprintf('Heat-trns coef:\tU = %0.2f\t\t[W/m2/K]\n',ST.U)
fprintf('Reqired area:\tA = %0.2f\t\t[m2]\n',ST.A)
fprintf('Adjusted area:\tAF = %0.2f\t\t[m2]\n',ST.AF)
fprintf('Pres. drop:\t\tdp = %0.2e\tdp = %0.2e\t[Pa]\n', H.ST.dp, C.ST.dp)
fprintf('Conv. coef.:\th = %0.2e\th = %0.2e\t[W/m2/K]\n', H.ST.h, C.ST.h)
fprintf('Nusselt number:\tNu = %0.1f\t\tNu = %0.1f\t\t[-]\n', H.ST.Nu, C.ST.Nu)
fprintf('HX Area:\t\tA = %0.2f\t\tA = %0.2f\t\t[m2]\n', H.ST.Aht, C.ST.Aht)
fprintf('Velocity:\t\tv = %0.2f\t\tv = %0.2f\t\t[m/s]\n\n', H.ST.v, C.ST.v)

fprintf('COUNTERFLOW PLATE HEAT EXCHANGER PROPERTIES\n')
fprintf('Heat-trns coef:\tU = %0.2f\t\t[W/m2/K]\n',CF.U)
fprintf('Reqired area:\tA = %0.2f\t\t[m2]\n',CF.A)
fprintf('Total length:\tL = %0.2f\t\t[m]\n',CF.L)
fprintf('Total width:\tW = %0.2f\t\t[m]\n',CF.W)
fprintf('Total height:\tH = %0.2f\t\t[m]\n',CF.H)
fprintf('Pres. drop:\t\tdp = %0.2e\tdp = %0.2e\t[Pa]\n', H.CF.dp, C.CF.dp)
fprintf('Conv. coef.:\th = %0.2e\th = %0.2e\t[W/m2/K]\n', H.CF.h, C.CF.h)
fprintf('Nusselt number:\tNu = %0.1f\t\tNu = %0.1f\t\t[-]\n', H.CF.Nu, C.CF.Nu)
fprintf('HX Area:\t\tA = %0.2f\t\tA = %0.2f\t\t[m2]\n', H.CF.Aht, C.CF.Aht)
fprintf('Velocity:\t\tv = %0.2f\t\tv = %0.2f\t\t[m/s]\n\n', H.CF.v, C.CF.v)

fprintf('COUNTERFLOW DUCTS HEAT EXCHANGER PROPERTIES\n')
fprintf('Heat-trns coef:\tU = %0.2f\t\t[W/m2/K]\n',CD.U)
fprintf('Reqired area:\tA = %0.2f\t\t[m2]\n',CD.A)
fprintf('Total length:\tL = %0.2f\t\t[m]\n',CD.L)
fprintf('Total width:\tW = %0.2f\t\t[m]\n',CD.W)
fprintf('Total height:\tH = %0.2f\t\t[m]\n',CD.Height)
fprintf('Pres. drop:\t\tdp = %0.2e\tdp = %0.2e\t[Pa]\n', CD.H.dp, CD.C.dp)
fprintf('Conv. coef.:\th = %0.2e\th = %0.2e\t[W/m2/K]\n', CD.H.h, CD.C.h)
fprintf('Nusselt number:\tNu = %0.1f\t\tNu = %0.1f\t\t[-]\n', CD.H.Nu, CD.C.Nu)
fprintf('HX Area:\t\tA = %0.2f\t\tA = %0.2f\t\t[m2]\n', CD.H.Aht, CD.C.Aht)
fprintf('Velocity:\t\tv = %0.2f\t\tv = %0.2f\t\t[m/s]\n', CD.H.v, CD.C.v)
fprintf('Reynolds:\t\tRe = %0.3e\tRe = %0.3e\t[-]\n\n', CD.H.Re, CD.C.Re)


% figure(1)
% plot([0,1],[C.Tout, C.Tin], [0,1],[H.Tin,H.Tout]);
% title('Counterflow configuration')
% grid on;
% legend('Cold', 'Hot');
% xlabel('Length [-]')
% ylabel('Temperature [K]')
