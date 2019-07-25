function [cc,p] = ccof(a,b)
[cc, p] = corrcoef(a,b,'rows','complete');
        if isnan(cc(2))
            disp('nan ccof');
        end
cc = cc(2); p = p(2);
end