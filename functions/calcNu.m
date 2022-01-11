function Nu = calcNu(dat)
n = 0.35;

if dat.Re > 1e4 && dat.Pr < 160 && dat.Pr > 0.7
    Nu = 0.023*dat.Re^(4/5)*dat.Pr^n;
elseif dat.Re < 2300 && dat.Pr < 16700 && dat.Pr > 0.48
    Nu = 1.86*(dat.De*dat.Re*dat.Pr/dat.L)^(1/3);
else
    Nu1 = 1.86*(dat.De*2300*dat.Pr/dat.L)^(1/3);
    Nu2 = 0.023*1e4^(4/5)*dat.Pr^n;
    Nu = Nu2 + (Nu1 - Nu2)*(1e4 - dat.Re)/(1e4 - 2300);
end

end