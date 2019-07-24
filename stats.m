function stats( i )
    if ~iscolumn(i); i = i'; end
    a.n = len(i);
    a.mean = mean(i);
    a.median=median(i);
    a.min=min(i);
    a.max=max(i);
    a.std=std(i);
    stat=a;
    stat
    %datastats(i)
end

