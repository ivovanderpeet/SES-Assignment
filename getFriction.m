function f = getFriction(dat)
    if dat.Re < 2300
        f = 64/dat.Re;
    elseif dat.Re > 2300 && dat.Re < 4000
        f = nan;
    else
        f = 0.0055*(1+(2e4*dat.rough/dat.De +1e6/dat.Re)^(1/3));
    end

end