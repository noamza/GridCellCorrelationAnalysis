function   shufflingMain()
    warning off;
    %%%
    binsize = 45;
    %fn = sprintf('C:\\Noam\\Data\\muscimol\\noam\\cells_%dmin_c_midall_gridscore.mat',binsize); first set
    fn = sprintf('C:\\Noam\\Data\\muscimol\\noam\\cells_%dmin_d_patchtraj_rayleigh',45);
    load(fn);
    [groups ~] = findSimultaneouslyRecordedCells(cells);
    gridThreshBef = 0.3; gridThreshMid = 0.25;    
    
    iCbmIuRbm250u = {}; 
    for ri = 1:length(groups)
        ri
        g = groups{ri}; %AM I ITERATIONG THROUGH ALL?
        for j = 1:length(g)-1
            for k = j+1:length(g)
                cb = shuffleTimeCorrelations (g(j).before,g(k).before, 1, 200, 0.006);
                cm = shuffleTimeCorrelations (g(j).midall,g(k).midall, 1, 200, 0.006);
                iu = Intersection_Over_Union(g(j).before.module.x, g(j).before.module.y, g(k).before.module.x,g(k).before.module.y);
                iCbmIuRbm250u{g(j).ind,g(k).ind} = [g(j).ind g(k).ind cb cm iu... 
                    g(j).before.rayleigh_score g(k).before.rayleigh_score g(j).midall.rayleigh_score g(k).midall.rayleigh_score ...
                    p];
            end
        end
    end
    
     iud = []; iun = []; iuo = []; ium = []; iua = []; iug = []; iudd = [];                                       
     for ri = 1:length(groups)
        g = groups{ri}; %NEED GROUP FOR NON DECREASING
        for j = 1:length(g)-1
            for k = j+1:length(g)
                i1 = g(j).ind;  i2 = g(k).ind; gm1 =  g(j).midall.gridscore; gm2 =  g(k).midall.gridscore;
                c1 = g(j).before; c2 = g(k).before; cb = iCbmIuRbm250u{i1,i2}(3); cm = iCbmIuRbm250u{i1,i2}(4);
                iu = Intersection_Over_Union(c1.module.x,c1.module.y, c2.module.x,c2.module.y);
                %t = [iu i1 i2 cb cm];
                t = iCbmIuRbm250u{i1,i2};
                %both cells are grid
                if c1.gridscore >= gridThreshBef  && c2.gridscore >= gridThreshBef
                        iug(end +1,:) = t;
                        %both decreasing
                    if     gm1 <= gridThreshMid  && gm2 <= gridThreshMid
                        iud(end +1,:) = t;
                        %both non decreasing
                    elseif gm1 >  gridThreshMid && gm2 >  gridThreshMid
                        iun(end +1,:) = t;
                        %one decreasing
                    elseif(gm1 <= gridThreshMid && gm2 >  gridThreshMid) ||...
                          (gm1 >  gridThreshMid && gm2 <  gridThreshMid)
                        iuo(end +1,:) = t;
                        %should be empty
                    else
                        ium(end +1,:) = t;
                    end
                end
                iua(end +1,:) = t;
            end
        end
     end
     iudd = [iun; iuo]; %grid but not decreasing
    %iCbmIuRbm250u              i1 i2 cb cm iu r1b r2b r1m r2m 
    %                           1  2  3  4  5  6   7   8   9
    i1i = 1; i2i = 2; cbi = 3; cmi = 4; iui = 5; r1bi = 6; r2bi = 7; r1mi = 8; r2mi = 9;
    
    % RAYLEIGH + MODULE
    criuBef = {}; criuMid = {};
    criuB =[]; criuM = [];
    for ri = 1:length(gids)
        gis = gids(ri)
        cis = cids{gis}
        g = groups{gis};
        g = g(cis);
        for j = 1:length(g)-1
            for k = j+1:length(g)
                c1 = g(j).before; c2 = g(k).before;
                c = shuffleTimeCorrelations (c1,c2, 1, 200, 0.006); 
                r1 = c1.rayleigh_score; r2 = c2.rayleigh_score;
                iu = Intersection_Over_Union(c1.module.x,c1.module.y, c2.module.x,c2.module.y);
                criuBef{g(j).ind,g(k).ind} = [c r1 r2 iu];
                criuB(end+1,:) = [g(j).ind g(k).ind c r1 r2 iu];
                c1 = g(j).midall; c2 = g(k).midall;
                c = shuffleTimeCorrelations (c1,c2, 1, 200, 0.006);
                r1 = c1.rayleigh_score; r2 = c2.rayleigh_score;
                iu = Intersection_Over_Union(c1.module.x,c1.module.y, c2.module.x,c2.module.y);
                criuMid{g(j).ind,g(k).ind} = [c r1 r2 iu];
                criuM(end+1,:) = [g(j).ind g(k).ind c r1 r2 iu];
            end
        end
    end

    %g(j).ind <= size(criuBef,1) && g(k).ind <= size(criuBef,2) && ~isempty(criuBef{g(j).ind,g(k).ind})
                    %t = criuBef{g(j).ind,g(k).ind};
    
    %Can you do three histograms of i/u, of %gridness-lost pairs 
                                            %gridness-not-lost pairs 
                                            %gridness-lost-and-unlost pairs?
     %CHECK 0's and NANS                                       

    %RAYLEIGH
    pd.gridThreshBef = gridThreshBef; pd.gridThreshMid = gridThreshMid;
    h=figure; si = 0; row = 4; 
    pd.sesh = 'MID'; pd.co = 'ro'; pd.string = 'decreasing gs pairs'; [h,si] = plotRow(iud,pd,h,1,si);
    %pd.sesh = 'MID'; pd.co = 'ro'; pd.string = 'non decreasing gs pairs';[h,si] = plotRow(iudd,pd,h,row,si);
    %pd.sesh = 'BEF'; pd.co = 'bo'; pd.string = 'decreasing gs pairs'; [h,si] = plotRow(iud,pd,h,row,si);
    %pd.sesh = 'BEF'; pd.co = 'bo'; pd.string = 'non decreasing gs pairs';[h,si] = plotRow(iudd,pd,h,row,si);
    
    
    naxes = 6; ncol = 2;
    nRow = ceil( naxes / ncol ) ;
    rowH = 0.75 / nRow ;  colW = 0.75 / ncol ; %width of plot
           %offset
    colX = 0.05 + linspace( 0, 0.9, ncol+1 ) ; colX = colX(1:end-1) ;
    rowY = 0.05 + linspace( 0.9, 0, nRow+1 ) ; rowY = rowY(2:end);    
    %rowId = ceil( ai / ncol );
    %colId = ai - (rowId - 1) * ncol;
    %x = colX(colId); y = rowY(rowId);
    %p = [x, y, colW, rowH];
     x = colX( ai - (ceil( ai / ncol ) - 1) * ncol); y = rowY(ceil( ai / ncol ));p = [x, y, colW, rowH];
    
    
    %RAYLEIGH FINAL
    row = 4; col = 2; ii =1; figure; a = [];
    t = iud; lims = [min(min(t(:,cbi), min(t(:,cmi)))) max(max(t(:,cbi)), max(t(:,cmi)))] ;
    tg = t(t(:,r2mi)>=0.5 & t(:,r1mi)>=0.5,:); tl = t(t(:,r2mi)<0.5 | t(:,r1mi)<0.5,:);  
    %r1 vs r2
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(t(:,r1mi),t(:,r2mi),'go'); 
    xlim([0 1]); ylim([0 1]); title(sprintf('decreasing gs pairs (%.2f>%.2f)\nrayleigh i1 vs i2 MID',gridThreshBef,gridThreshMid));  
    ii = ii+1;
    %cb vs cr less than cluster
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(tl(:,cbi),tl(:,cmi),'ro'); title('rayleigh < 0.5 corr b vs m');
    hold on; plot([min(tl(:,cbi)) max(tl(:,cbi))], polyval( polyfit( tl(:,cbi),tl(:,cmi),1),[min(tl(:,cbi)) max(tl(:,cbi))] ) ); 
    xlim(lims); ylim(lims);
    %cb vs cr greater than cluster
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(tg(:,cbi),tg(:,cmi),'ro'); title('rayleigh > 0.5 corr b vs m');
    hold on; plot([min(tg(:,cbi)) max(tg(:,cbi))], polyval( polyfit( tg(:,cbi),tg(:,cmi),1),[min(tg(:,cbi)) max(tg(:,cbi))] ) );
    xlim(lims); ylim(lims);
    %hists
    bn = -0.05:0.005:0.05; 
    a(end+1) = subplot(row,col,ii);ii=ii+1;hist(tl(:,cbi),bn);title(sprintf('%s\n%s','bef: hist corrs','~(r1&r2 > 0.5)'));xlim([-0.05 0.05]);
    a(end+1) = subplot(row,col,ii);ii=ii+1;hist(tg(:,cbi),bn);title(sprintf('%s\n%s','bef: hist corrs',' (r1&r2 > 0.5)'));xlim([-0.05 0.05]);
    a(end+1) = subplot(row,col,ii);ii=ii+1;hist(tl(:,cmi),bn);title(sprintf('%s\n%s','mid: hist corrs','~(r1&r2 > 0.5)'));xlim([-0.05 0.05]);
    a(end+1) = subplot(row,col,ii);ii=ii+1;hist(tg(:,cmi),bn);title(sprintf('%s\n%s','mid: hist corrs',' (r1&r2 > 0.5)'));xlim([-0.05 0.05]);
    %end
    
    
    f = figure(2); tp = uix.HBox( 'Parent', f, 'Spacing', 0);lp = uix.VBox('Parent', tp );rp = uix.VBox('Parent', tp);
    set( tp, 'Widths', [-1,-1] );% Adjust the main layout
    %vc= uicontainer('Parent', lp);% t = [t -1];
     %RAYLEIGH FINAL
    row = 4; col = 2; ii =1; a = [];
    t = iud; lims = [min(min(t(:,cbi), min(t(:,cmi)))) max(max(t(:,cbi)), max(t(:,cmi)))] ;
    tg = t(t(:,r2mi)>=0.5 & t(:,r1mi)>=0.5,:); tl = t(t(:,r2mi)<0.5 | t(:,r1mi)<0.5,:);  
    %r1 vs r2
    a(end+1) = axes('Parent',uicontainer('parent',lp)); plot(a(end),t(:,r1mi),t(:,r2mi),'go'); 
    xlim(a(end),[0 1]); ylim(a(end),[0 1]); title(a(end),sprintf('decreasing gs pairs (%.2f>%.2f)\nrayleigh i1 vs i2 MID',gridThreshBef,gridThreshMid));
    %cb vs cr less than cluster
    a(end+1) = axes('Parent',uicontainer('parent',lp)); plot(a(end),tl(:,cbi),tl(:,cmi),'ro'); title(a(end),'rayleigh < 0.5 corr b vs m');
    hold(a(end),'on'); plot(a(end),[min(tl(:,cbi)) max(tl(:,cbi))], polyval( polyfit( tl(:,cbi),tl(:,cmi),1),[min(tl(:,cbi)) max(tl(:,cbi))] ) ); 
    xlim(a(end),lims); ylim(a(end),lims);
    %cb vs cr greater than cluster
    a(end+1) = axes('Parent',uicontainer('parent',lp)); plot(a(end),tg(:,cbi),tg(:,cmi),'ro'); title(a(end),'rayleigh > 0.5 corr b vs m');
    hold(a(end),'on'); plot(a(end),[min(tg(:,cbi)) max(tg(:,cbi))], polyval( polyfit( tg(:,cbi),tg(:,cmi),1),[min(tg(:,cbi)) max(tg(:,cbi))] ) );
    xlim(a(end),lims); ylim(a(end),lims);
    %hists
    bn = -0.05:0.005:0.05;
    rp1 = uix.VBox('Parent', rp,'Padding', 0, 'Spacing', 0 );%title(a(end),sprintf('%s\n%s','bef: hist corrs','~(r1&r2 > 0.5)'));title(a(end),sprintf('%s\n%s','bef: hist corrs',' (r1&r2 > 0.5)'));
    a(end+1) = axes('Parent',uicontainer('parent',rp1)); hist(a(end),tl(:,cbi),bn);xlim(a(end),[-0.05 0.05]); set(a(end),'xticklabel',[]);%xlabel(a(end),[]); 
    a(end+1) = axes('Parent',uicontainer('parent',rp1)); hist(a(end),tg(:,cbi),bn);xlim(a(end),[-0.05 0.05]);
    rp2 = uix.VBox('Parent', rp,'Padding', 0, 'Spacing', 0 );%title(a(end),sprintf('%s\n%s','mid: hist corrs','~(r1&r2 > 0.5)'));title(a(end),sprintf('%s\n%s','mid: hist corrs',' (r1&r2 > 0.5)'));
    a(end+1) = axes('Parent',uicontainer('parent',rp2)); hist(a(end),tl(:,cmi),bn);xlim(a(end),[-0.05 0.05]);set(a(end),'xticklabel',[]);%axis(a(end),'xticklabel',[]); 
    a(end+1) = axes('Parent',uicontainer('parent',rp2)); hist(a(end),tg(:,cmi),bn);xlim(a(end),[-0.05 0.05]);
    set( rp, 'heights', [-1,-1] );
    %end
    
    
    
    
    
    viewLayout = uix.VBox( 'Parent', lp,'Padding', 0, 'Spacing', 0);
    t = [];
    set(viewLayout, 'Heights', t);
    
    
    
    a(end+1) = subplot(row,col,ii);ii=ii+1; hist(t(:,cbi)); title('hist corrs bef');
    a(end+1) = subplot(row,col,ii);ii=ii+1; hist(t(:,cmi)); title('hist corrs mid'); 

    
    
    
    figure; %hist iu
    subplot(1,3,1);hist(iud,20); title('IU pairs decreasing'); xlabel('intersection / union');
    subplot(1,3,2);hist(iun(iun>0),20); title('IU pairs non-decreasing');
    subplot(1,3,3);hist(iua(iua>0),20); title('IU pairs all'); 
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SHUFFLING
    [gids, cids] = groupsByGridThresh(groups, gridThreshBef, gridThreshMid); % >grid, <mid
    gcorrBefore = {}; gcorrMid = {};
     corrBefore = {};  corrMid = {};
    n = 250;
    tic
    for ri = 1:length(gids)
        gis = gids(ri)
        cis = cids{gis}
        g = groups{gis};
        g = g(cis);
        for j = 1:length(g)-1
            for k = j+1:length(g)
                [ri j k] 
                gcorrBefore{k,j,ri} = shuffleTimeCorrelations (g(j).before, g(k).before, n, 200, 0.006);
                gcorrMid{k,j,ri}    = shuffleTimeCorrelations (g(j).midall, g(k).midall, n, 200, 0.006);
                corrBefore{g(j).ind,g(k).ind} = gcorrBefore{k,j,ri};
                   corrMid{g(j).ind,g(k).ind} = gcorrMid{k,j,ri};
            end
        end
    end
    disp('************************');
    toc
      
    
%     t = corrBefore{11,12};
%     [s,I] = sort(abs(t),'descend');
%     p = round(find(I==1)/length(I),2) %index of non shuffled 
    ri = 11, ci = 12, c1 = g(1).before, c2 = g(2).before, 
    lag = 200, sigma = 0.006;
    
    ri = 11, ci = 12
    [r, c] = size(corrBefore);
    pvalsBefMid = cell(r,c); pbs = []; pms = []; ijcbmpbm = []; 
    for ri = 1:r
        for ci = 1:c
            if ~isempty(corrBefore{ri,ci})
                t = corrBefore{ri,ci}; cb = t(1);
                [s,I] = sort(abs(t),'descend');
                pb = round(find(I==1)/length(I),2); %index of non shuffled
                t = corrMid{ri,ci}; cm = t(1);
                [s,I] = sort(abs(t),'descend');
                pm = round(find(I==1)/length(I),2);
                pvalsBefMid{ri,ci} = [pb pm]; pbs(end+1) = pb; pms(end+1) = pm;
                ijcbmpbm(end+1,:) = [ri ci cb cm pb pm ];
            end
        end
    end
    figure(1);
    subplot(2,1,1);hist(pbs*100,100); title('p shuffle before 250 update');
    subplot(2,1,2);hist(pms*100,100); title('p shuffle mid 250 update');
    
    i1i = 1; i2i = 2; cbi = 3; cmi = 4; pbi = 5; pmi = 6; 
    
    %HIST C B M
    figure; bn = [-.05:0.0025:.05]; %LEARN DIFFERENCES BETWEEN HIST AND HISTOGRAM
    subplot(2,1,1);histogram(ijcbmpbm(:,cbi),bn); title('corr by shuffle(250) before'); hold on
    histogram(ijcbmpbm(ijcbmpbm(:,pbi)<=0.01,cbi),bn,'FaceColor','r','EdgeColor','r');
    legend({'corrs for pairs above gs 0.3 bef 0.25 mid','pval above 0.01'});
    subplot(2,1,2);histogram(ijcbmpbm(:,cmi),bn); title('corr by shuffle(250) middle'); hold on
    histogram(ijcbmpbm(ijcbmpbm(:,pmi)<=0.01,cmi),bn,'FaceColor','r','EdgeColor','r');
    %SPARSE
    
    %CORR CORR P
    figure; plot(ijcbmpbm(:,cbi),ijcbmpbm(:,cmi),'.'); hold on;
    plot(ijcbmpbm(ijcbmpbm(:,pmi)<=0.01,cbi),ijcbmpbm(ijcbmpbm(:,pmi)<=0.01,cmi),'o');
    %plot([min(ijcbmpbm(:,cbi)) max(ijcbmpbm(:,cbi))], polyval( polyfit( ijcbmpbm(:,cbi),ijcbmpbm(:,cmi),1),[min(ijcbmpbm(:,cbi)) max(ijcbmpbm(:,cbi))]));
    f1=fit(ijcbmpbm(:,cbi),ijcbmpbm(:,cmi),'poly1');plot(f1,'b');
    legend({'corr b vs m','significant mid corr <= 0.01',sprintf('%.2fx + %.2f',round(f1.p1,2),round(f1.p2,2))});
    axis equal; ylim([-.05 .15]); xlim([-.05 .15]); xlabel(''); ylabel('')
    title('corr before vs middle by shuffling p');
    
    
    figure
    f1=fit(ijcbmpbm(:,cbi),ijcbmpbm(:,cmi),'poly1');%fit(cdate,pop,'poly1')
    plot(f1);
    legend(sprintf('%.2fx+%.2f',round(f1.p1,2),round(f1.p2,2)));
    
    disp('');
    %}
end





