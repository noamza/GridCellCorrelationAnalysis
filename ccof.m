function [cc,p] = ccof(a,b)
[cc, p] = corrcoef(a,b);
        if isnan(cc(2))
            disp('nan ccof');
        end
cc = cc(2); p = p(2);
end