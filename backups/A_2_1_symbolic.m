clear all; close all; clc;
%% Constants
THin = 240;     % [C] Given
TCout = 110;    % [C] Given
Qdot = 5e7;     % [W] Given
mdotC = 150;    % [kg/s] Given

%% Loop for determining CpC

% Initial Guess
CpC = 4200;     % [J/kg/K] Assumed constant

% Iteratively determine CpC
ERR = 1;
while ERR > 1e-10
    CC = mdotC*CpC;
    TCin = TCout - Qdot/(mdotC*CpC); % [K] Calculate cold-side inlet temperature
    TmC = (TCin + TCout)/2;
    CpC_ = XSteam('CpL_T', TmC)*1000;
    ERR = abs(CpC_ - CpC)/CpC;

    CpC = CpC_;
end

% Final value for TCin and TmC
TCin = TCout - Qdot/(mdotC*CpC); % [K] Calculate cold-side inlet temperature
TmC = (TCin + TCout)/2;

%% Symbolic hot-side stuff
syms CH

CR = CH/CC; % [-] Known that CH < Cc
THout = THin - CC./CH.*(TCout-TCin);

% For counterflow
dT1 = THout - TCin;
dT2 = THin - TCout;
dTlm = (dT2 - dT1)./log(dT2./dT1);

UA = Qdot./dTlm;
NTU = UA./CH; % [-] Known that CH < Cc

% Determine unknowns from NTU == 3
CH =    double(vpasolve(NTU==3, CH));
CR =    double(subs(CR));
THout = double(subs(THout));
dT1 =   double(subs(dT1));
dTlm =  double(subs(dTlm));
UA =    double(subs(UA));
NTU =   double(subs(NTU));

% Determine CpH from THout and THin
TmH = (THout + THin)/2;
CpH = XSteam('CpL_T', TmH)*1000;

% Calcualte mdotH from CH and CpH
mdotH = CH/CpH;

% Calculate effectiveness
eps = (1-exp(-NTU.*(1-CR)))./(1-CR.*exp(-NTU.*(1-CR)));

% Saturation pressure at given temperature
pH = XSteam('psat_T', TmH);
pC = XSteam('psat_T', TmC);

% Determine mean density
rhoH = XSteam('rhoL_T', TmH);
rhoC = XSteam('rhoL_T', TmC);

% Determine average dynamic viscosity (assumed saturated liquid p & T)
muH = XSteam('my_pT', pH, TmH*1.001); % T multiplied by 1.001 to ensure liquid domain
muC = XSteam('my_pT', pC, TmC*1.001);

% Determine thermal conductivity
kH = XSteam('tcL_T', TmH);
kC = XSteam('tcL_T', TmC);

%% Figures & display
fprintf('Heat capacity:\tCpH = %0.2f,\t\tCpC = %0.2f\t[J/kg/K]\n', CpH, CpC)
fprintf('Mass flows:\t\tmdotH = %3.2f,\t\tmdotC = %3.2f\t[kg/s]\n', mdotH, mdotC)
fprintf('Density:\t\trhoH = %0.2f,\t\tfrhoC = %0.2f\t[kg/m3]\n',rhoH, rhoC)
fprintf('Dyn. visc.:\t\tmuH = %0.3e,\tmuC = %0.3e\t[Pa s]\n',muH, muC)
fprintf('Nondimensional:\tNTU = %0.2f,\t\t\teps = %0.4f\n', NTU, eps)
fprintf('Heat transfer:\tUA = %0.3e\t[W/m2/K]\n',UA)

figure(1)
plot([0,1],[TCout, TCin], [0,1],[THin,THout]);
title('Counterflow configuration')
grid on;
legend('Cold', 'Hot');
xlabel('Length [-]')
ylabel('Temperature [K]')
