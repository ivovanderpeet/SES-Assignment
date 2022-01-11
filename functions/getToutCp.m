function dat = getToutCp(dat)
bool_mdot = ~max(contains(fieldnames(dat),'mdot')) && max(contains(fieldnames(dat),'V'));
bool_Tin = max(contains(fieldnames(dat),'Tin'));
bool_Tout = max(contains(fieldnames(dat),'Tout'));

dat.Cp = 4200; % initial guess;
dat.rho = 1000; % initial guess;
if bool_mdot
    dat.mdot = dat.V*dat.rho;
end


ERR = 1;
while ERR > 1e-10
    if bool_Tin
        dat.Tout = dat.Qdot/dat.mdot/dat.Cp + dat.Tin;
    elseif bool_Tout
        dat.Tin = -dat.Qdot/dat.mdot/dat.Cp + dat.Tout;
    end
    dat.Tm = (dat.Tout + dat.Tin)/2;
    if bool_mdot
        dat.rho = XSteam('rho_pT', dat.p/1e5, dat.Tm);
        dat.mdot = dat.V*dat.rho;
    end
    Cp_ = XSteam('Cp_pT', dat.p/1e5, dat.Tm)*1e3;
    ERR = abs(Cp_ - dat.Cp)/dat.Cp;
    dat.Cp = Cp_;
end

dat.Tout = dat.Qdot/dat.mdot/dat.Cp + dat.Tin;
dat.Tm = (dat.Tout + dat.Tin)/2;
dat.C = dat.mdot*dat.Cp;
end