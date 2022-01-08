function [H,C] = counterPlate(H,C)
%% Hot side
% Flow properties
H.CF.A = H.CF.H*H.CF.W;
H.CF.Atot = H.CF.A*H.CF.Nplate; %% OF MOET DIT 1200 ZIJN?
H.CF.v = H.V/H.CF.Atot;

PH = 2*H.CF.H + 2*H.CF.W;
H.CF.De = 4*H.CF.A/PH;
H.CF.Re = H.rho*H.CF.v*H.CF.De/H.mu;

% Heat transfer properties
H.CF.Aht = (H.CF.Nplate*2-1)*H.CF.W*H.CF.L;
H.CF.Pr = H.Pr;
H.CF.Nu = calcNu(H.CF);
H.CF.h = H.CF.Nu*H.k/H.CF.De;

% Pressure drop
H.CF.f = getFriction(H.CF);
H.CF.dp = H.CF.f*H.CF.L/H.CF.De *H.rho*H.CF.v^2/2;

%% Cold side
% Flow properties
C.CF.A = C.CF.H*C.CF.W;
C.CF.Atot = C.CF.A*C.CF.Nplate;
C.CF.v = C.V/C.CF.Atot;

PC = 2*C.CF.H + 2*C.CF.W;
C.CF.De = 4*C.CF.A/PC;
C.CF.Re = C.rho*C.CF.v*C.CF.De/C.mu;

% Heat transfer properties
C.CF.Aht = (C.CF.Nplate*2-1)*C.CF.W*C.CF.L;
C.CF.Pr = C.Pr;
C.CF.Nu = calcNu(C.CF);
C.CF.h = C.CF.Nu*C.k/C.CF.De;

% Pressure drop
C.CF.f = getFriction(C.CF);
C.CF.dp = C.CF.f*C.CF.L/C.CF.De *C.rho*C.CF.v^2/2;

end