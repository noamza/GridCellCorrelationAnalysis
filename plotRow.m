function [h,si] = plotRow(data,pd,h,r,si)
    d = plotRowShuffRayData(data,pd);
    [h,si] = plotRowShuffRay(d,h,r,length(d),si);
    
end

function [h,si] = plotRowShuffRay(pdata,h,r,c,si)
    i = 0; ax = []; h = figure(h);
    
    i=i+1; d = pdata{i}; si=si+1;
    ax(end+1) = subplot(r,c,si); hold on; 
    plot(ax(end),d.x, d.y, d.co); 
    xlim(ax(end),d.xl); ylim(ax(end),d.yl); title(ax(end),d.title);
    
    i=i+1; d = pdata{i}; si=si+1;
    ax(end+1) = subplot(r,c,si); hold on; 
    plot(ax(end),d.x, d.y, d.co); title(ax(end),d.title);
    plot(ax(end),[min(d.x) max(d.x)], polyval( polyfit( d.x, d.y,1), [min(d.x) max(d.x)] ) );
    xlim(ax(end),d.xl); ylim(ax(end),d.yl);
    
    i=i+1; d = pdata{i}; si=si+1;
    ax(end+1) = subplot(r,c,si); hold on; 
    plot(ax(end),d.x, d.y, d.co); title(ax(end),d.title);
    plot(ax(end),[min(d.x) max(d.x)], polyval( polyfit( d.x, d.y, 1), [min(d.x) max(d.x)] ) );
    xlim(ax(end),d.xl); ylim(ax(end),d.yl);
    
    i=i+1; d = pdata{i}; si=si+1;
    ax(end+1) = subplot(r,c,si); hold on; 
    hold on; hist(ax(end),d.x); 
    xlim(ax(end),d.xl); title(ax(end),d.title);
    
    i=i+1; d = pdata{i}; si=si+1;
    ax(end+1) = subplot(r,c,si); hold on; 
    hold on; hist(ax(end),d.x); 
    xlim(ax(end),d.xl); title(ax(end),d.title);
    
    i=i+1; d = pdata{i}; si=si+1;
    ax(end+1) = subplot(r,c,si); hold on; 
    hold on; hist(ax(end),d.x); 
    xlim(ax(end),d.xl); title(ax(end),d.title);
    
    i=i+1; d = pdata{i}; si=si+1;
    ax(end+1) = subplot(r,c,si); hold on; 
    hold on; hist(ax(end),d.x); 
    xlim(ax(end),d.xl); title(ax(end),d.title);
end

function pdata = plotRowShuffRayData(data,pd)
    i1i = 1; i2i = 2; cbi = 3; cmi = 4; iui = 5; r1bi = 6; r2bi = 7; r1mi = 8; r2mi = 9;
    t = data;
    
    if isequal(pd.sesh,'MID')
        r2 = r2mi;
        r1 = r1mi;
    else
        r2 = r2bi;
        r1 = r1bi;
    end
    
    tg = t(t(:,r2)>=0.5 & t(:,r1)>=0.5,:); tl = t(t(:,r2)<0.5 | t(:,r1)<0.5,:);
    lims = [min(min(t(:,cbi), min(t(:,cmi)))) max(max(t(:,cbi)), max(t(:,cmi)))] ;
    
    ri = 0; pdr = {}; %MID BEF
    
    %1) r1 vs r2 clustering
    r=[];r.title=sprintf('[%s] %s (%.2f>%.2f)\nrayleigh i1 vs i2',pd.sesh, pd.string,pd.gridThreshBef,pd.gridThreshMid);
    r.x = t(:,r1); r.y = t(:,r2); r.co = pd.co; r.xl=[0 1]; r.yl=[0 1]; 
    ri=ri+1; pdr{ri}=r;
    %2) corrs of lower cluster
    r=[];r.title=sprintf('%s','rayleigh < 0.5 corr b vs m');
    r.x = tl(:,cbi); r.y = tl(:,cmi); r.co = pd.co; r.xl=lims; r.yl=lims; 
    ri=ri+1; pdr{ri}=r;
    %3) corrs of higher cluster
    r=[];r.title=sprintf('%s','rayleigh > 0.5 corr b vs m');
    r.x = tg(:,cbi); r.y = tg(:,cmi); r.co = pd.co; r.xl=lims; r.yl=lims; 
    ri=ri+1; pdr{ri}=r;
    
    %4) 5) do hists before mid by cluster (x4) 
    
    r=[];r.title=sprintf('%s\n%s','bef: hist corrs','~(r1&r2 > 0.5)');
    r.x = tl(:,cbi);r.xl=[-0.05 0.15]; 
    ri=ri+1; pdr{ri}=r;
    
    r=[];r.title=sprintf('%s\n%s','bef: hist corrs','r1&r2 > 0.5');
    r.x = tg(:,cbi);r.xl=[-0.05 0.15]; 
    ri=ri+1; pdr{ri}=r;
    
    r=[];r.title=sprintf('%s\n%s','mid: hist corrs','~(r1&r2 > 0.5)');
    r.x = tl(:,cmi);r.xl=[-0.05 0.15]; 
    ri=ri+1; pdr{ri}=r;
    pdata = pdr;
    
    r=[];r.title=sprintf('%s\n%s','mid: hist corrs','r1&r2 > 0.5');
    r.x = tg(:,cmi);r.xl=[-0.05 0.15]; 
    ri=ri+1; pdr{ri}=r;
    pdata = pdr;
end
%{
    t = iud; lims = [min(min(t(:,cbi), min(t(:,cmi)))) max(max(t(:,cbi)), max(t(:,cmi)))] ;
    tg = t(t(:,r2mi)>=0.5 & t(:,r1mi)>=0.5,:); tl = t(t(:,r2mi)<0.5 | t(:,r1mi)<0.5,:);
    
a(end+1) = subplot(row,col,ii);ii=ii+1; plot(t(:,r1mi),t(:,r2mi),'ro'); 
    xlim([0 1]); ylim([0 1]); title(sprintf('decreasing gs pairs (%.2f>%.2f)\nrayleigh i1 vs i2 MID',gridThreshBef,gridThreshMid));  
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(tg(:,cmi),tg(:,cbi),'ro'); title('rayleigh > 0.5 corr b vs m');
    hold on; plot([min(tg(:,cbi)) max(tg(:,cbi))], polyval( polyfit( tg(:,cmi),tg(:,cbi),1),[min(tg(:,cbi)) max(tg(:,cbi))] ) );
    xlim(lims); ylim(lims);
    
a(end+1) = subplot(row,col,ii);ii=ii+1; plot(tl(:,cmi),tl(:,cbi),'ro'); title('rayleigh < 0.5 corr b vs m');
    hold on; plot([min(tl(:,cbi)) max(tl(:,cbi))], polyval( polyfit( tl(:,cmi),tl(:,cbi),1),[min(tl(:,cbi)) max(tl(:,cbi))] ) ); 
    xlim(lims); ylim(lims);
    
a(end+1) = subplot(row,col,ii);ii=ii+1; hist(t(:,cbi)); title('hist corrs bef');
    
a(end+1) = subplot(row,col,ii);ii=ii+1; hist(t(:,cmi)); title('hist corrs mid');



    %INTERSECTION UNION
    row = 5; col = 3; ii = 1; 
    figure; 
     t = iud; t = t(t(:,iui)>0,:); tg = t(t(:,iui)>=0.6,:); tl = t(t(:,iui)<0.6,:); 
    a(end+1) = subplot(row,col,ii);ii=ii+1; hist(t(:,iui),20); 
    title(sprintf('decreasing gs pairs (%.2f>%.2f) i/u',gridThreshBef,gridThreshMid)); xlim([0 1]);
    a(end+1) = subplot(row,col,ii);ii=ii+1;plot(tg(:,cmi),tg(:,cbi),'ro');title('iu >= 0.6 corr b vs m');
    hold on; plot([min(tg(:,cbi)) max(tg(:,cbi))], polyval( polyfit( tg(:,cmi),tg(:,cbi),1),[min(tg(:,cbi)) max(tg(:,cbi))] ) );
    axis equal;
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(tl(:,cmi),tl(:,cbi),'bo'); title('iu < 0.6 corr b vs m');
    hold on; plot([min(tl(:,cbi)) max(tl(:,cbi))], polyval( polyfit( tl(:,cmi),tl(:,cbi),1),[min(tl(:,cbi)) max(tl(:,cbi))] ) ); 
    axis equal;
    t = iun;t = t(t(:,iui)>0,:); tg = t(t(:,iui)>=0.6,:); tl = t(t(:,iui)<0.6,:); 
    a(end+1) = subplot(row,col,ii);ii=ii+1; hist(t(:,iui),20); title('non decreasing pairs i/u');xlim([0 1]);
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(tg(:,cmi),tg(:,cbi),'ro'); title('iu >= 0.6 corr b vs m');
    hold on; plot([min(tg(:,cbi)) max(tg(:,cbi))], polyval( polyfit( tg(:,cmi),tg(:,cbi),1),[min(tg(:,cbi)) max(tg(:,cbi))] ) );
    axis equal;
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(tl(:,cmi),tl(:,cbi),'bo'); title('iu < 0.6 corr b vs m');
    hold on; plot([min(tl(:,cbi)) max(tl(:,cbi))], polyval( polyfit( tl(:,cmi),tl(:,cbi),1),[min(tl(:,cbi)) max(tl(:,cbi))] ) );  
    axis equal;
     t = iuo; t = t(t(:,iui)>0,:); tg = t(t(:,iui)>=0.6,:); tl = t(t(:,iui)<0.6,:); 
    a(end+1) = subplot(row,col,ii);ii=ii+1; hist(t(:,iui),20); title('one decreasing pairs i/u');xlim([0 1]);
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(tg(:,cmi),tg(:,cbi),'ro'); title('iu >= 0.6 corr b vs m');
    hold on; plot([min(tg(:,cbi)) max(tg(:,cbi))], polyval( polyfit( tg(:,cmi),tg(:,cbi),1),[min(tg(:,cbi)) max(tg(:,cbi))] ) );
    axis equal;
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(tl(:,cmi),tl(:,cbi),'bo'); title('iu < 0.6 corr b vs m');
    hold on; plot([min(tl(:,cbi)) max(tl(:,cbi))], polyval( polyfit( tl(:,cmi),tl(:,cbi),1),[min(tl(:,cbi)) max(tl(:,cbi))] ) );
    axis equal;
     t = ium; t = t(t(:,iui)>0,:); tg = t(t(:,iui)>=0.6,:); tl = t(t(:,iui)<0.6,:); 
    a(end+1) = subplot(row,col,ii);ii=ii+1; hist(t(:,iui),20); title('rest of pairs i/u');xlim([0 1]);
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(tg(:,cmi),tg(:,cbi),'ro'); title('iu >= 0.6 corr b vs m');
    hold on; plot([min(tg(:,cbi)) max(tg(:,cbi))], polyval( polyfit( tg(:,cmi),tg(:,cbi),1),[min(tg(:,cbi)) max(tg(:,cbi))] ) );   
    axis equal;
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(tl(:,cmi),tl(:,cbi),'bo'); title('iu < 0.6 corr b vs m');
    hold on; plot([min(tl(:,cbi)) max(tl(:,cbi))], polyval( polyfit( tl(:,cmi),tl(:,cbi),1),[min(tl(:,cbi)) max(tl(:,cbi))] ) );
    axis equal;
     t = iua; t = t(t(:,iui)>0,:); tg = t(t(:,iui)>=0.6,:); tl = t(t(:,iui)<0.6,:); 
    a(end+1) = subplot(row,col,ii);ii=ii+1; hist(t(:,iui),20); title('all pairs i/u');xlim([0 1]);
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(tg(:,cmi),tg(:,cbi),'ro'); title('iu >= 0.6 corr b vs m');
    hold on; plot([min(tg(:,cbi)) max(tg(:,cbi))], polyval( polyfit( tg(:,cmi),tg(:,cbi),1),[min(tg(:,cbi)) max(tg(:,cbi))] ) );
    axis equal;
    a(end+1) = subplot(row,col,ii);         plot(tl(:,cmi),tl(:,cbi),'bo'); title('iu < 0.6 corr b vs m');
    hold on; plot([min(tl(:,cbi)) max(tl(:,cbi))], polyval( polyfit( tl(:,cmi),tl(:,cbi),1),[min(tl(:,cbi)) max(tl(:,cbi))] ) );
    axis equal;


    
    t = findObj(gcf,'type','axes');
    for ri = 1:length(t)
        a = t(ri);
        if mod(ri-1,3) == 0
            mod(ri-1,3)
            axis(a, 'equal');
        end
    end
%}