function ax = plotTR(c,str,ax)%SETS LIMITS OF ARENA TO 100
    title(ax, str.tt);xlabel(ax, str.xt);ylabel(ax, str.yt);hold(ax, 'on');
    
    plot(ax, c.px, c.py);
    scatter(ax, c.sx, c.sy, '.');
    xlim(ax, [0 100]), ylim(ax, [0 100]);
    
    axis(ax,'tight'); set(ax,'ydir','normal');colormap(ax,'jet');
end
