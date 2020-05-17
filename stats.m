function stat = stats( a,nans,infs)
    a = toCol(a);
    if nargin > 1 
        a=a(~isnan(a)); % only non-nans
        if nargin==3; a=a(a~=inf&a~=-inf);end %only non infs
    end
    s.n = len(a);
    s.mean = mean(a);
    s.median=median(a);
    s.min=min(a);
    s.max=max(a);
    s.std=std(a);
    stat=s;
    %datastats(i)
end

