function ax = plotCCSpace(c1,c2,str,ax,p)%SETS LIMITS OF ARENA TO 100 
    title(ax, str.tt);xlabel(ax, str.xt);ylabel(ax, str.yt);hold(ax, 'on');
    
    cc = xcorr2(c2.rm,c1.rm); %reverse order for perspective
    imagesc(ax, cc); 
    plot(ax, size(cc,2)/2,size(cc,1)/2,'md','MarkerFaceColor','w','MarkerSize',7);
    %intersection union
    if isfield(p.module)
    plot(ax, c1.module.x, c1.module.y,'w'); plot(ax, c2.module.x, c2.module.y,'r');
    title(ax, sprintf('Int/Union=%.2f',Intersection_Over_Union(...
        c2.module.x,c2.module.y,c1.module.x,c1.module.y) ),'fontsize',9);
    end
    
   axis(ax,'tight'); set(ax,'ydir','normal');colormap(ax,'jet');
end
