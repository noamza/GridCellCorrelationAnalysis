function ax = plotTR(ax,c,str,p)%SETS LIMITS OF ARENA TO 100
    title(ax, str.t);
    xlabel(ax, str.x);
    ylabel(ax, str.y);
    hold(ax, 'on');
    plot(ax, c.px, c.py,'k','linewidth',1.5);
    hold(ax, 'on');
    plot(ax, c.sx, c.sy, 'r.','markersize',2.7),...
%     plot(ax, c.px, c.py);
%     scatter(ax, c.sx, c.sy);
    %xlim(ax, [0 100]), ylim(ax, [0 100]);
    axis(ax,'tight'); 
    set(ax,'ydir','normal');colormap(ax,'jet');
    axis(ax,'square'); 
    if exist('p','var') && (isfield(p,'on'))
        %set(ax, 'XTick', []); set(ax, 'YTick', []);
    else
        axis(ax,'off');
        set(ax, 'XTick', []); set(ax, 'YTick', []);
    end
    %set(ax,p);
end

