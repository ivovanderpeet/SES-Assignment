function [H,C] = shellTube(H,C)
%% Hot side
% Flow properties
H.ST.A = pi/4*H.ST.ID^2;
H.ST.Atot = H.ST.A*H.ST.Ntube; %% OF MOET DIT 1200 ZIJN?
H.ST.v = H.V/H.ST.Atot;

PHI = pi*H.ST.ID; % Inner perimeter
PHO = pi*H.ST.OD; % Outer preimeter
H.ST.De = 4*H.ST.A/PHI;
H.ST.Re = H.rho*H.ST.v*H.ST.ID/H.mu;

% Heat transfer properties
H.ST.Aht = PHI*H.ST.L*H.ST.Ntube*H.ST.Npass;
H.ST.Pr = H.Pr;
H.ST.Nu = calcNu(H.ST);
H.ST.h = H.ST.Nu*H.k/H.ST.De;

% Pressure drop
H.ST.f = getFriction(H.ST);
H.ST.dp = H.ST.f*H.ST.L*H.ST.Npass/H.ST.De *H.rho*H.ST.v^2/2;

%% Cold side
% Flow properties
C.ST.A = pi/4*C.ST.ID^2 - H.ST.Atot*H.ST.Npass;
C.ST.v = C.V/C.ST.A;

PCht = PHO*H.ST.Ntube*H.ST.Npass; % Heat exchanging perimeter
PC = pi*C.ST.ID + PCht;        % Total perimeter
C.ST.De = 4*C.ST.A/PC;
C.ST.Re = C.rho*C.ST.v*C.ST.De/C.mu;

% Heat transfer properties
C.ST.Aht = PCht*C.ST.L;
C.ST.Pr = C.Pr;
C.ST.Nu = calcNu(C.ST);
C.ST.h = C.ST.Nu*C.k/C.ST.De;

% Pressure drop
C.ST.f = getFriction(C.ST);
C.ST.dp = C.ST.f*C.ST.L/C.ST.De *C.rho*C.ST.v^2/2;

end