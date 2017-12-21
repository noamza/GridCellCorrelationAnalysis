function ax = plotHD(c,str,ax)%SETS LIMITS OF ARENA TO 100
    title(ax, str.tt);xlabel(ax, str.xt);ylabel(ax, str.yt);hold(ax, 'on');
    
    bins= 0:3:360;
    c.si = discretize(c.st, [-Inf; mean([c.pt(2:end) c.pt(1:end-1)],2); +Inf]);
    rd = histcounts(rad2deg(c.hd(c.si))+180,bins)./(histcounts(rad2deg(c.hd)+180,bins)+0.001);
    rd = smooth(rd,15);
    plot(ax,bins(2:end),rd);

    axis(ax,'tight'); set(ax,'ydir','normal');colormap(ax,'jet');
end
