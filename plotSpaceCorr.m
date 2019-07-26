function ax = plotSpaceCorr(ax,c1,c2,p,str)%SETS LIMITS OF ARENA TO 100 
    title(ax, str.t);xlabel(ax, str.x);
    ylabel(ax, str.y);hold(ax, 'on');
    
    %cc = xcorr2(c2.rm,c1.rm); %reverse order for perspective
    cc = xcorr2g(c2.rm,c1.rm); %reverse order for perspective
    %cc = xcorr2(c2.rm-mean(c2.rm(:)),c1.rm-mean(c1.rm(:))); %reverse order for perspective
    cc = imgaussfilt(cc,1.4);
    cenr = round(size(cc,2)/2); cenc = round(size(cc,1)/2);
    imagesc(ax, cc);
    %text(cenr,cenc,sprintf('%.1f',cc(cenr,cenc)),...'Color','w');%'FontSize',10);
    %plot(ax, size(cc,2)/2,size(cc,1)/2,'md','MarkerFaceColor','w','MarkerSize',7); %mark 0,0
    %intersection union
    if isfield(p,'module')
        plot(ax, c1.module.x, c1.module.y,'w'); 
        plot(ax, c2.module.x, c2.module.y,'r');
        title(ax, sprintf('Int/Union=%.2f',Intersection_Over_Union(...
        c2.module.x,c2.module.y,c1.module.x,c1.module.y) ),'fontsize',9);
    end
        
   axis(ax,'off'); axis(ax,'square');
   axis(ax,'tight'); set(ax,'ydir','normal');
   colormap(ax,'jet');
    if (isfield(p,'off'))
        set(ax, 'XTick', []); set(ax, 'YTick', []);
    end
end
