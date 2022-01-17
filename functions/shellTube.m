function [T,S] = shellTube(T,S)
%% Hot side
% Flow properties
T.A = pi/4*T.ID^2;
T.Atot = T.A*T.Ntube; %% OF MOET DIT 1200 ZIJN?
T.v = T.V/T.Atot;

PHI = pi*T.ID; % Inner perimeter
PHO = pi*T.OD; % Outer preimeter
T.De = 4*T.A/PHI;
T.Re = T.rho*T.v*T.ID/T.mu;

% Heat transfer properties
T.Aht = PHI*T.L*T.Ntube*T.Npass;
T.Pr = T.Pr;
T.Nu = calcNu(T);
T.h = T.Nu*T.k/T.De;

% Pressure drop
T.f = getFriction(T);
T.dp = T.f*T.L*T.Npass/T.De *T.rho*T.v^2/2;

%% Cold side
% Flow properties
S.A = pi/4*S.ID^2 - T.Atot*T.Npass;
S.v = S.V/S.A;

PCht = PHO*T.Ntube*T.Npass; % Heat exchanging perimeter
PC = pi*S.ID + PCht;        % Total perimeter
S.De = 4*S.A/PC;
S.Re = S.rho*S.v*S.De/S.mu;

% Heat transfer properties
S.Aht = PCht*S.L;
S.Pr = S.Pr;
S.Nu = calcNu(S);
S.h = S.Nu*S.k/S.De;

% Pressure drop
S.f = getFriction(S);
S.dp = S.f*S.L/S.De *S.rho*S.v^2/2;

end