function h =  plotTimeSpaceCC(p)
    binl = p.binspike;%0.06 0.1;%.3
    lag = p.lag;
    nbins = round(360/p.bindeg); %rounds good!  fprintf('nbins %d\n', nbins);
    sigma = p.sigma;
    parent = p.parent;
    h = p.fig;
    
end