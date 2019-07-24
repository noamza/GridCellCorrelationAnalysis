function ax = plotHD(ax,c,str)%SETS LIMITS OF ARENA TO 100
    title(ax, str.t);xlabel(ax, str.x);ylabel(ax, str.y);hold(ax, 'on');
    bins= 0:.1:360;
    c.si = discretize(c.st, [-Inf; mean([c.pt(2:end) c.pt(1:end-1)],2); +Inf]);
    rd = histcounts(wrapTo360(rad2deg(c.hd(c.si))),bins)./...
        (histcounts(wrapTo360(rad2deg(c.hd)),bins)+0.00001);
    rd = movmean(rd,round(len(rd)/5));
    p=plot(ax,bins(2:end),rd/1); p.LineWidth = 1.5; 
    set(ax,'xtick',[0, 360],'xticklabel',{sprintf('0%c',char(176));...
       sprintf('360%c',char(176))},'yticklabel',{},'ylim',[0 max(rd)]) 
    axis(ax, 'square');colormap(ax,'jet');
    %figure, plot(rd/sum(rd)), ylim([0 1/10000])
end