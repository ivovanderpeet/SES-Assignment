function sat = getSatCurve()
T = linspace(-10,374,1000);
len = length(T);
sL = zeros(size(T));
sV = sL;
hV = sL;
hL = sL;
rhoL = sL;
rhoV = sL;
p = sL;
for i = 1:len
    sL(i) = XSteam('sL_T', T(i))*1e3;
    sV(i) = XSteam('sV_T', T(i))*1e3;
    hL(i) = XSteam('hL_T', T(i))*1e3;
    hV(i) = XSteam('hV_T', T(i))*1e3;
    p(i) = XSteam('psat_T',T(i))*1e5;
    rhoL(i) = XSteam('rhoL_T',T(i));
    rhoV(i) = XSteam('rhoV_T',T(i));
end

sat.p = [p(~isnan(hL)),flip(p(~isnan(hV)))];
sat.T = [T(~isnan(hL)),flip(T(~isnan(hV)))];
sat.s = [sL(~isnan(sL)),flip(sV(~isnan(sV)))];
sat.h = [hL(~isnan(hL)),flip(hV(~isnan(hV)))];
sat.v = [1./rhoL(~isnan(hL)),flip(1./rhoV(~isnan(hV)))];

end