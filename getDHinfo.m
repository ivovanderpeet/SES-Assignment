function dat = getDHinfo(dat)
% Initial Guess
dat.Cp = 4200;     % [J/kg/K] Assumed constant
dat.mdot = dat.V*1000;

% Iteratively determine dat.Cp
ERR = 1;
while ERR > 1e-10
    dat.C = dat.mdot * dat.Cp;
    dat.Tin = dat.Tout - dat.Qdot/(dat.mdot*dat.Cp); % [K] Calculate cold-side inlet temperature
    dat.Tm = (dat.Tin + dat.Tout)/2;
    Cp_ = XSteam('Cp_pT', dat.p/1e5, dat.Tm)*1000;
    dat.rho = XSteam('rho_pT', dat.p/1e5, dat.Tm);
    dat.mdot = dat.V*dat.rho;
    ERR = abs(Cp_ - dat.Cp)/dat.Cp;

    dat.Cp = Cp_;
end
clear('ERR');

% Final value for dat.Tin and dat.Tm
dat.Tin = dat.Tout - dat.Qdot/(dat.mdot*dat.Cp); % [K] Calculate cold-side inlet temperature
dat.Tm = (dat.Tin + dat.Tout)/2;
dat.mdot = 150*XSteam('rho_pT', dat.p/1e5, dat.Tm)/1000;
end