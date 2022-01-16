function [H,C] = HX(H,C)
%% Hot side
% Flow properties
Rcond.HX.A = Rcond.HX.D^2;
Rcond.HX.Atot = Rcond.HX.A*(Rcond.HX.Nplatex*Rcond.HX.Nplatey*2); %% OF MOET DIT 1200 ZIJN?
Rcond.HX.v = Rcond.V/Rcond.HX.Atot;

PH = 4*Rcond.HX.D; % Inner perimeter
Rcond.HX.De = 4*Rcond.HX.A/PH;
Rcond.HX.Re = Rcond.rho*Rcond.HX.v*Rcond.HX.De/Rcond.mu;

% Heat transfer properties
Rcond.HX.Aht = PH*Rcond.HX.L*(Rcond.HX.Nplatex*Rcond.HX.Nplatey*2 - (Rcond.HX.Nplatex + Rcond.HX.Nplatey)/2);%(Rcond.HX.Nplatex*2-1+Rcond.HX.Nplatey*2-1)*Rcond.HX.W*Rcond.HX.L;
Rcond.HX.Pr = Rcond.Pr;
Rcond.HX.Nu = calcNu(Rcond.HX);
Rcond.HX.h = Rcond.HX.Nu*Rcond.k/Rcond.HX.De;

% Pressure drop
Rcond.HX.f = getFriction(Rcond.HX);
Rcond.HX.dp = Rcond.HX.f*Rcond.HX.L/Rcond.HX.De*Rcond.rho*Rcond.HX.v^2/2;

%% Cold side
% Flow properties
% Flow properties
DHcond.HX.A = DHcond.HX.D^2;
DHcond.HX.Atot = DHcond.HX.A*(DHcond.HX.Nplatex*DHcond.HX.Nplatey*2); %% OF MOET DIT 1200 ZIJN?
DHcond.HX.v = DHcond.V/DHcond.HX.Atot;

PC = 4*DHcond.HX.D; % Inner perimeter
DHcond.HX.De = 4*DHcond.HX.A/PC;
DHcond.HX.Re = DHcond.rho*DHcond.HX.v*DHcond.HX.De/DHcond.mu;

% Heat transfer properties
DHcond.HX.Aht = DHcond.HX.L*PC*(DHcond.HX.Nplatex*DHcond.HX.Nplatey*2 - (DHcond.HX.Nplatex+DHcond.HX.Nplatey)/2);%(DHcond.HX.Nplatex*2-1+DHcond.HX.Nplatey*2-1)*DHcond.HX.W*DHcond.HX.L;
DHcond.HX.Pr = DHcond.Pr;
DHcond.HX.Nu = calcNu(DHcond.HX);
DHcond.HX.h = DHcond.HX.Nu*DHcond.k/DHcond.HX.De;

% Pressure drop
DHcond.HX.f = getFriction(DHcond.HX);
DHcond.HX.dp = DHcond.HX.f*DHcond.HX.L/DHcond.HX.De*DHcond.rho*DHcond.HX.v^2/2;

end