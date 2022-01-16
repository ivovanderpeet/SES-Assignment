function [H,C] = counterDucts(H,C)
%% Hot side
% Flow properties
H.A = H.D^2;
H.Atot = H.A*(H.Nplatex*H.Nplatey*2);
H.v = H.V/H.Atot;

PH = 4*H.D;
H.De = 4*H.A/PH;
H.Re = H.rho*H.v*H.De/H.mu;

% Heat transfer properties
H.Aht = PH*H.L*(H.Nplatex*H.Nplatey*2 - (H.Nplatex + H.Nplatey)/2);
H.Pr = H.Pr;
H.Nu = calcNu(H);
H.h = H.Nu*H.k/H.De;

% Pressure drop
H.f = getFriction(H);
H.dp = H.f*H.L/H.De *H.rho*H.v^2/2;

%% Cold side
% Flow properties
C.A = C.D^2;
C.Atot = C.A*(C.Nplatex*C.Nplatey*2);
C.v = C.V/C.Atot;

PC = 4*C.D;
C.De = 4*C.A/PC;
C.Re = C.rho*C.v*C.De/C.mu;

% Heat transfer properties
C.Aht = C.L*PC*(C.Nplatex*C.Nplatey*2 - (C.Nplatex+C.Nplatey)/2);
C.Pr = C.Pr;
C.Nu = calcNu(C);
C.h = C.Nu*C.k/C.De;

% Pressure drop
C.f = getFriction(C);
C.dp = C.f*C.L/C.De *C.rho*C.v^2/2;

end