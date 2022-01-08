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

%% Heat exchanger
H.ID = 20e-3;
H.OD = H.ID + 2* 2e-3;
H.pitch = 1.25*H.OD; % Distance between centers of two tubes 1.25 OD is conventional
H.Ntube = 210;
H.Npass = 2;
H.L = 5;
H.rough = 0.1e-3;
C.rough = H.rough;

% https://www.engineeringtoolbox.com/smaller-circles-in-larger-circle-d_1849.html
% http://hydra.nat.uni-magdeburg.de/packing/cci/
load('circ.mat');
C.ID = H.pitch*circ.ratio(H.Ntube*H.Npass);
C.L = H.L;

% Unsure about C.h (effective diameter)
[H,C] = shellTube(H,C);

% Correction factor Shell & Tube
Z = (H.Tin-H.Tout)/(C.Tout-C.Tin);
Y = (H.Tout-H.Tin)/(C.Tin-H.Tin);
F = 0.93;

kWall = 50; % Steel = 50 W/m/K
U = (1/C.h + H.OD*log(H.OD/H.ID)/2/kWall + H.OD/H.ID/H.h)^-1;  % 4PC00 7.43
A = UA/U;
AF = A/F;



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
fprintf('Velocity:\t\tv = %0.2f\t\tv = %0.2f\t\t[m/s]\n\n', H.v, C.v)

fprintf('OVERALL & NONDIMENSIONAL PROPERTIES\n')
fprintf('Heat transfer:\tUA = %0.3e\t[W/K]\n',UA)
fprintf('Transfer units:\tNTU = %0.2f\t\t[-]\n',NTU)
fprintf('Effectiveness:\teps = %0.4f\t[-]\n', eps)
fprintf('Effectiveness2:\teps2 = %0.4f\t[-]\n\n', eps2)

fprintf('HEAT EXCHANGER PROPERTIES\n')
fprintf('Corr. factor:\tZ = %0.3f\t\t[-]\n',Z)
fprintf('Corr. factor:\tY = %0.3f\t\t[-]\n',Y)
fprintf('Corr. factor:\tF = %0.3f\t\t[-]\n',F)
fprintf('Heat-trns coef:\tU = %0.2f\t\t[W/m2/K]\n',U)
fprintf('Reqired area:\tA = %0.2f\t\t[m2]\n',A)
fprintf('Adjusted area:\tAF = %0.2f\t\t[m2]\n',AF)
fprintf('Outer area:\t\tAO = %0.2f\t\t[m2]\n',C.Aht)
fprintf('Inner area:\t\tAI = %0.2f\t\t[m2]\n',H.Aht)


% figure(1)
% plot([0,1],[C.Tout, C.Tin], [0,1],[H.Tin,H.Tout]);
% title('Counterflow configuration')
% grid on;
% legend('Cold', 'Hot');
% xlabel('Length [-]')
% ylabel('Temperature [K]')
