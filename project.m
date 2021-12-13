%% Constants

T1 = 240;
T2 = 110;
mdot2 = 0.15; % m3/s
cp = 4200;
Qdot = 50*10^6;

%% A.2.1 Design of a heat-only plant

% 1. Calculate the return temperature of the district-heating stream via eq.(2.6).

T3 = (mdot2*cp*T2 - Qdot)/mdot2*cp

%%
