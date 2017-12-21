%c session
%str strings for labels
%ax axis
function ax = plotRM(c,str,ax)%m,n,l    
    title(ax, str.tt);xlabel(ax, str.xt);ylabel(ax, str.yt);hold(ax, 'on');
    
    imagesc(ax, imgaussfilt(c.rm,2,'FilterDomain','spatial')); 
    
    axis(ax,'tight'); set(ax,'ydir','normal');colormap(ax,'jet'); axis(ax,'square');
end

