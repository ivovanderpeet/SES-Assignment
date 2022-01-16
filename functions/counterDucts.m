function [x] = counterDucts(x)
%% Hot side
% Flow properties
x.A = x.D^2;
x.Atot = x.A*(x.Nplatex*x.Nplatey*2);
x.v = x.V/x.Atot;

PH = 4*x.D;
x.De = 4*x.A/PH;
x.Re = x.rho*x.v*x.De/x.mu;

% Heat transfer properties
x.Aht = PH*x.L*(x.Nplatex*x.Nplatey*2 - (x.Nplatex + x.Nplatey)/2);
x.Pr = x.Pr;
x.Nu = calcNu(x);
x.h = x.Nu*x.k/x.De;

% Pressure drop
x.f = getFriction(x);
x.dp = x.f*x.L/x.De *x.rho*x.v^2/2;
end