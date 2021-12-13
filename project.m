%% Constants

T_SGin = 240 + 273; % K
T_Cout = 110 + 273; % K
mdot = 150;     % kg/s
cp = 4200;      % J/K
Qdot = 50*10^6; % J/s

%% A.2.1 Design of a heat-only plant

% 1. Calculate the return temperature of the district-heating stream via eq.(2.6).

T_Cin = (mdot*cp*T_Cout - Qdot)/(mdot*cp)

% 3. Choose a value for the capacity flow ratio CR, see § 2.5 (rationale?)

CR = ;

% 4. Calculate the exit temperature of the well water and well mass flow using the relations in Sec. 2.5.
%    Judge whether the well mass flow is “reasonable” by comparison with the district-heating mass flow and adjust your choice for CR if deemed “unreasonable”.1



%%
