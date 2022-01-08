function [H,C] = counterDucts(H,C)
%% Hot side
% Flow properties
H.CD.A = H.CD.H*H.CD.W;
H.CD.Atot = H.CD.A*H.CD.Nplatex*H.CD.Nplatey; %% OF MOET DIT 1200 ZIJN?
H.CD.v = H.V/H.CD.Atot;

PH = 2*H.CD.H + 2*H.CD.W; % Inner perimeter
H.CD.De = 4*H.CD.A/PH;
H.CD.Re = H.rho*H.CD.v*H.CD.De/H.mu;

% Heat transfer properties
H.CD.Aht = (H.CD.Nplatex*2-1+H.CD.Nplatey*2-1)*H.CD.W*H.CD.L;
H.CD.Pr = H.Pr;
H.CD.Nu = calcNu(H.CD);
H.CD.h = H.CD.Nu*H.k/H.CD.De;

% Pressure drop
H.CD.f = getFriction(H.CD);
H.CD.dp = H.CD.f*H.CD.L/H.CD.De *H.rho*H.CD.v^2/2;

%% Cold side
% Flow properties
C.CD.A = C.CD.H*C.CD.W;
C.CD.Atot = C.CD.A*C.CD.Nplatex*C.CD.Nplatey;
C.CD.v = C.V/C.CD.Atot;

PC = 2*C.CD.H + 2*C.CD.W; % Inner perimeter
C.CD.De = 4*C.CD.A/PC;
C.CD.Re = C.rho*C.CD.v*C.CD.De/C.mu;

% Heat transfer properties
C.CD.Aht = (C.CD.Nplatex*2-1+C.CD.Nplatey*2-1)*C.CD.W*C.CD.L;
C.CD.Pr = C.Pr;
C.CD.Nu = calcNu(C.CD);
C.CD.h = C.CD.Nu*C.k/C.CD.De;

% Pressure drop
C.CD.f = getFriction(C.CD);
C.CD.dp = C.CD.f*C.CD.L/C.CD.De *C.rho*C.CD.v^2/2;

end