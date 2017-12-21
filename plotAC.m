function ax = plotAC(c,str,ax)
    title(ax, str.tt);xlabel(ax, str.xt);ylabel(ax, str.yt); hold(ax, 'on');
    
    imagesc(ax, imgaussfilt(c.ac,2,'FilterDomain','spatial'));
    x=xlim(ax); y=ylim(ax);
    %title(ax, sprintf('%sspk%d',s, length(p.sx)));
    lw = 0.5;
    %center point
    if length(c.module.hex_peaks) == 7 %CHANGE TO EXISTS
            plot(ax, c.module.x, c.module.y,'w','LineWidth',lw); 
            plot(ax, c.module.hex_peaks(:,1),c.module.hex_peaks(:,2),'wo','LineWidth',lw);
    end
    xlim(ax, x); ylim(ax, y);
    
    axis(ax,'tight'); set(ax,'ydir','normal');colormap(ax,'jet');
end

