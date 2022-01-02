function [H,C] = shellTube(H,C)
%% Hot side
% Flow properties
H.A = pi/4*H.ID^2;
H.Atot = H.A*H.Ntube;
H.V = H.mdot/H.rho;
H.v = H.V/H.Atot;

PHI = pi*H.ID; % Inner perimeter
PHO = pi*H.OD; % Outer preimeter
H.De = 4*H.A/PHI;
H.Re = H.rho*H.v*H.ID/H.mu;

% Heat transfer properties
H.Aht = PHI*H.L*H.Ntube*H.Npass;
H.Nu = calcNu(H);
H.h = H.Nu*H.k/H.De;

%% Cold side
% Flow properties
C.A = pi/4*C.ID^2 - H.Atot*H.Npass;
C.V = C.mdot/C.rho;
C.v = C.V/C.A;

PCht = PHO*H.Ntube*H.Npass; % Heat exchanging perimeter
PC = pi*C.ID + PCht;        % Total perimeter
C.De = 4*C.A/PC;
C.Re = C.rho*C.v*C.De/C.mu;

% Heat transfer properties
C.Aht = PCht*C.L;
C.Nu = calcNu(C);
C.h = C.Nu*C.k/C.De;

end