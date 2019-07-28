function [r,p,s] = ccof(a,b,rd)
s='';
[r, p] = corrcoef(a,b,'rows','complete');
if isnan(r(2))
    disp('nan ccof');
end
r = r(2); p = p(2);
if nargin==3
    f=['%0.' n2(rd) 'f'];
    s=sprintf(['r=' f ' p=' f],rnd(r,rd),rnd(p,rd));
end
end