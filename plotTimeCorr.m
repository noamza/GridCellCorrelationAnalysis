function ax = plotTimeCorr(c1,c2,params,str,ax)
    title(ax, str.tt);xlabel(ax, str.xt);ylabel(ax, str.yt); hold(ax, 'on');
    [p, ~] = compareByMovingDirection(c1,c2, params);
    X = p{1,2}'; Y = p{1,1}';
    plot(ax, X, Y,'linewidth', 3); hold(ax, 'on'); 
    plot(ax,  [0 0],  [-1 1],'r','linewidth', 1);
    text(ax, 0,0, sprintf('%.3f %.3f', round(Y(round(length(Y)/2)),3)), 'Color','r'); 
    axis(ax, 'tight');axis(ax, 'square');ax.Color = 'black'; colormap(ax,'jet');
    ylim(ax, [-max(abs(Y(:))) max(abs(Y(:)))]);
    ax.YAxis.Exponent = 0; ax.XAxis.TickLabelFormat = '%gs';
    
end

