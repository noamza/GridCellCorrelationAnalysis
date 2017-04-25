function y = round2r(x,r)
    y = floor(x) + round((x-floor(x))/r)*r;