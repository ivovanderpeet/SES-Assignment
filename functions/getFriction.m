function f = getFriction(dat)
    if dat.Re < 2300
        f = 64/dat.Re;
    elseif dat.Re > 2300 && dat.Re < 3000
        f = nan;
    else
        A = (2.457*log((7/dat.Re)^.9 + 0.27*dat.rough/dat.De))^16;
        B = (37530/dat.Re)^16;
        f = 8*((8/dat.Re)^12 + 1/(A+B)^1.5)^(1/12);
%         f = 0.0055*(1+(2e4*dat.rough/dat.De +1e6/dat.Re)^(1/3));
    end

end