function [datK,datU] = getToutCp(datK,datU)
datU.Cp = 4200; % initial guess;
ERR = 1;
while ERR > 1e-10
    datU.Tout = -datK.Qdot/datU.mdot/datU.Cp+ datU.Tin;
    datU.Tm = (datU.Tout + datU.Tin)/2;
    Cp_ = XSteam('Cp_pT', datU.p/1e5, datU.Tm)*1e3;
    ERR = abs(Cp_ - datU.Cp)/datU.Cp;
    datU.Cp = Cp_;
end

datU.Tout = -datK.Qdot/datU.mdot/datU.Cp + datU.Tin;
datU.Tm = (datU.Tout + datU.Tin)/2;
end