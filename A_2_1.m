clear all; close all; clc;
%% Constants
H.Tin = 240;     % [C] Given
C.Tout = 110;    % [C] Given
Qdot = 5e7;      % [W] Given
C.mdot = 150;    % [kg/s] Given

%% Loop for determining C.Cp

% Initial Guess
C.Cp = 4200;     % [J/kg/K] Assumed constant

% Iteratively determine C.Cp
ERR = 1;
while ERR > 1e-10
    C.C = C.mdot*C.Cp;
    C.Tin = C.Tout - Qdot/(C.mdot*C.Cp); % [K] Calculate cold-side inlet temperature
    C.Tm = (C.Tin + C.Tout)/2;
    C.Cp_ = XSteam('CpL_T', C.Tm)*1000;
    ERR = abs(C.Cp_ - C.Cp)/C.Cp;

    C.Cp = C.Cp_;
end
clear('ERR');

% Final value for C.Tin and C.Tm
C.Tin = C.Tout - Qdot/(C.mdot*C.Cp); % [K] Calculate cold-side inlet temperature
C.Tm = (C.Tin + C.Tout)/2;

%% Symbolic hot-side stuff
syms CH
H.C = CH;

CR = C.C/H.C; % [-] Given that C.C < H.C
H.Tout = H.Tin - C.C./H.C.*(C.Tout-C.Tin);

% For counterflow
dT1 = H.Tout - C.Tin;
dT2 = H.Tin - C.Tout;
dTlm = (dT2 - dT1)./log(dT2./dT1);

UA = Qdot./dTlm;
NTU = UA./C.C; % [-] Given that C.C < H.C

% Determine unknowns from CR = 1
H.C = double(vpasolve(CR==1, H.C));
CH = H.C;
CR =    double(subs(CR));
H.Tout = double(subs(H.Tout));
dT1 =   double(subs(dT1));

if CR == 1
    dTlm = dT1;
    UA = Qdot/dTlm;
    NTU = UA/C.C;
else
    dTlm =  double(subs(dTlm));
    UA =    double(subs(UA));
    NTU =   double(subs(NTU));
end

clear('CH');

% Determine H.Cp from H.Tout and H.Tin
H.Tm = (H.Tout + H.Tin)/2;
H.Cp = XSteam('CpL_T', H.Tm)*1000;

% Calcualte H.mdot from H.C and H.Cp
H.mdot = H.C/H.Cp;

% Calculate effectiveness
if CR == 1
    eps = NTU/(NTU+1);
else
    eps = (1-exp(-NTU.*(1-CR)))./(1-CR.*exp(-NTU.*(1-CR)));
end

% Saturation pressure at given temperature
H.psat = XSteam('psat_T', H.Tm);
C.psat = XSteam('psat_T', C.Tm);

% Determine mean density
H.rho = XSteam('rhoL_T', H.Tm);
C.rho = XSteam('rhoL_T', C.Tm);

% Determine average dynamic viscosity (assumed saturated liquid p & T)
H.mu = XSteam('my_pT', H.psat, H.Tm*.99999); % T multiplied by .99 to ensure liquid domain
C.mu = XSteam('my_pT', C.psat, C.Tm*.999);

% Determine thermal conductivity
H.k = XSteam('tcL_T', H.Tm);
C.k = XSteam('tcL_T', C.Tm);

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
H.ST.rough = 0.1e-3;
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
H.CF.Nplate = 25;
C.CF.Nplate = H.CF.Nplate;

H.CF.rough = 0.1e-3;
C.CF.rough = H.CF.rough;

[H,C] = counterPlate(H,C);

% Heat transfer coef.
CF.t = 2e-3;
CF.kwall = 50; % STEEL
CF.U = (1/H.CF.h + CF.t/CF.kwall + 1/C.CF.h)^-1;
CF.A = UA/CF.U;

%% Counterflow ducts heat exchanger
H.CD.L = 3;
C.CD.L = H.CD.L;
H.CD.H = 10e-3;
C.CD.H = H.CD.H;
H.CD.W = 10e-3;
C.CD.W = H.CD.W;
H.CD.Nplatey = 25;
C.CD.Nplatey = H.CD.Nplatey;
H.CD.Nplatex = 50;
C.CD.Nplatex = H.CD.Nplatex;

H.CD.rough = 0.1e-3;
C.CD.rough = H.CD.rough;

[H,C] = counterDucts(H,C);

% Heat transfer coef.
CD.t = 2e-3;
CD.kwall = 50; % STEEL
CD.U = (1/H.CD.h + CD.t/CD.kwall + 1/C.CD.h)^-1;
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
fprintf('Therm. cond.:\tk = %0.4f\t\tk = %0.4f\t\t[W/m/K]\n',H.k, C.k)
fprintf('Velocity:\t\tv = %0.2f\t\tv = %0.2f\t\t[m/s]\n\n', H.ST.v, C.ST.v)

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
fprintf('Outer area:\t\tAO = %0.2f\t\t[m2]\n',C.ST.Aht)
fprintf('Inner area:\t\tAI = %0.2f\t\t[m2]\n\n',H.ST.Aht)

fprintf('COUNTERFLOW PLATE HEAT EXCHANGER PROPERTIES\n')
fprintf('Heat-trns coef:\tU = %0.2f\t\t[W/m2/K]\n',CF.U)
fprintf('Reqired area:\tA = %0.2f\t\t[m2]\n',CF.A)
fprintf('Outer area:\t\tAO = %0.2f\t\t[m2]\n',C.CF.Aht)
fprintf('Inner area:\t\tAI = %0.2f\t\t[m2]\n',H.CF.Aht)

fprintf('COUNTERFLOW DUCTS HEAT EXCHANGER PROPERTIES\n')
fprintf('Heat-trns coef:\tU = %0.2f\t\t[W/m2/K]\n',CD.U)
fprintf('Reqired area:\tA = %0.2f\t\t[m2]\n',CD.A)
fprintf('Outer area:\t\tAO = %0.2f\t\t[m2]\n',C.CD.Aht)
fprintf('Inner area:\t\tAI = %0.2f\t\t[m2]\n',H.CD.Aht)

% figure(1)
% plot([0,1],[C.Tout, C.Tin], [0,1],[H.Tin,H.Tout]);
% title('Counterflow configuration')
% grid on;
% legend('Cold', 'Hot');
% xlabel('Length [-]')
% ylabel('Temperature [K]')
