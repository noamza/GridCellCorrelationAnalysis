warning off;
%%%
binsize = 45;
%fn = sprintf('C:\\Noam\\Data\\muscimol\\noam\\cells_%dmin_c_midall_gridscore.mat',binsize); first set
fn = sprintf('C:\\Noam\\Data\\muscimol\\noam\\cells_%dmin_d_patchtraj_rayleigh',45);
load(fn);
[groups ~] = findSimultaneouslyRecordedCells(cells);
gridThreshBef = 0.3; gridThreshMid = 0.25;        
[gids, cids] = groupsByGridThresh(groups, gridThreshBef, gridThreshMid);
%iterate through all groups
%correlate all pairs
%smooth at different sized windows
ijcb = []; ijcm = [];%corrBefore = {};  corrMid = {}; 
tic
for ri = 1:length(gids)
    ri
    %gis = gids(ri); cis = cids{gids(ri)};%gis}
    g = groups{gids(ri)};%gis};
    g = g(cids{gids(ri)});%cis);
    for j = 1:length(g)-1
        for k = j+1:length(g)
            %BEFORE
            bcorsmooth = []; mcorsmooth = [];
            c1 = g(j).before; c2 = g(k).before;
            spkt1=floor(c1.st*1000)+1; spkt2=floor(c2.st*1000)+1; min_time=min(min(spkt1),min(spkt2)); %TIME IN MS
            offset = 100; spkt1 = spkt1-min_time+offset; spkt2 = spkt2-min_time+offset; max_time=max(max(spkt1),max(spkt2))+offset;
            train1=zeros(1,max_time);train2 =zeros(1,max_time); train1(spkt1)=1; train2(spkt2)=1; %LOTS OF 0'ss
            for e = 0:16            
                train1smooth = movmean(train1,2^e);
                train2smooth = movmean(train2,2^e);
                bcorsmooth(end+1) = corr(train1smooth', train2smooth'); %xcov(train1smooth, train2smooth,0, 'coef');
            end
            %DURING
            c1 = g(j).midall; c2 = g(k).midall;
            spkt1=floor(c1.st*1000)+1; spkt2=floor(c2.st*1000)+1; min_time=min(min(spkt1),min(spkt2)); %TIME IN MS
            offset = 100; spkt1 = spkt1-min_time+offset; spkt2 = spkt2-min_time+offset; max_time=max(max(spkt1),max(spkt2))+offset;
            train1=zeros(1,max_time);train2 =zeros(1,max_time); train1(spkt1)=1; train2(spkt2)=1; %LO
            for e = 0:16            
                train1smooth = movmean(train1,2^e);
                train2smooth = movmean(train2,2^e);
                mcorsmooth(end+1) = corr(train1smooth', train2smooth');%xcov(train1smooth, train2smooth,0, 'coef');
            end
            %gcorrBefore{k,j,ri} = shuffleTimeCorrelations (g(j).before, g(k).before, n, 200, 0.006);
            %gcorrMid{k,j,ri}     = shuffleTimeCorrelations (g(j).midall, g(k).midall, n, 200, 0.006);
            %corrBefore{g(j).ind,g(k).ind} = bcorsmooth;
            %corrMid{g(j).ind,g(k).ind}     = mcorsmooth;
            ijcb(end+1,:)  = [g(j).ind,g(k).ind, bcorsmooth];
            ijcm(end+1,: )= [g(j).ind,g(k).ind, mcorsmooth];
        end
    end
end
toc


figure;
for i = 1:17
    subplot(3,6,i);
    x = ijcb(:,i+2); y = ijcm(:,i+2); a = 0; b = 1.6;
    y(x>b*std(x)) = a; x(x>b*std(x)) = a; x(y>b*std(y)) = a; y(y>b*std(y)) = a; x = x(x~=0); y = y(y~=0); %go by std of x
    plot(x,y,'.'); axis('tight');a = xlim; b = ylim; xlim([min(a(1),b(1)), max(a(2),b(2))]); ylim(xlim)   %'mo','markersize',16,'linewidth',lfs
    hold on;
    title(sprintf('%dms',2^(i-1)));
    p = polyfit(x,y,1);
    yfit = polyval(p,x);yfit =  p(1) * x + p(2);
    plot(x,yfit)
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq = 1 - SSresid/SStotal;
    %legend((sprintf('r = %.1f %.1f', polyfit( x, y,1) )  ),'Location','best');
    legend( sprintf('r = %.1f res = %.2f',p(1),rsq)  ,'Location','best');
    hold off;%axis('equal');%axis('square');
end
suptitle('std 1.6');

    %plot([min(x) max(x)], polyval( polyfit( x, y,1), [min(x) max(x)] )); %'linewidth',lfs,'color','b'
    %b = x\y;
    %plot(x,b*x);
    %text(gca,0.8*max(d.x),polyval(polyfit(d.x, d.y, 1),0.8*max(d.x)), sprintf('%.2f %.2f', round(polyfit(d.x, d.y, 1),2)),'Color','b','FontSize',afs);
    %text(sprintf('%.1f %.1f', round( polyfit(x, y, 1), 1 ) ) );   
    
%do spatial smoothing, sample middle do before and after

ijspacb = []; ijspacm = []; e = 6.^(-1:0.2:1); e = [1e-9 e];%corrBefore = {};  corrMid = {}; 
tic
%figure; i = 0;
for ri = 1:length(gids)
    ri
    %gis = gids(ri); cis = cids{gids(ri)};%gis}
    g = groups{gids(ri)};%gis};
    g = g(cids{gids(ri)});%cis);
    for j = 1:length(g)-1
        for k = j+1:length(g)
            tb = []; tm = [];
            for i = 1:length(e);
                sigma = e(i);
                %i = i+1; subplot(12,10,i);
                c1 = g(j).before; c2 = g(k).before;
                cc = xcorr2(c2.rm,c1.rm); %reverse order for perspective
                cenr = round(size(cc,2)/2); cenc = round(size(cc,1)/2);
                ccg = imgaussfilt(cc,sigma); imagesc(ccg); title(ccg(cenr,cenc));
                tb(end+1) = [ccg(cenr,cenc)]; 
                %imagesc(cc); %xlabel(az, sprintf('c%d x c%d',c1.ind,c2.ind),'fontweight','bold');colormap(az,'jet');
                %imagesc(imgaussfilt(cc,1)); title(ccg(cenr,cenc))%colormap('jet');
                %MIDALL
                c1 = g(j).midall; c2 = g(k).midall;
                cc = xcorr2(c2.rm,c1.rm); %reverse order for perspective
                cenr = round(size(cc,2)/2); cenc = round(size(cc,1)/2);
                ccg = imgaussfilt(cc,sigma);
                tm(end+1) = [ccg(cenr,cenc)];
            end
            ijspacb(end+1,:)  = [g(j).ind,g(k).ind, tb];
            ijspacm(end+1,: )= [g(j).ind,g(k).ind, tm];
        end
    end
end
toc


figure;
for i = 1:12
    subplot(3,4,i);
    sigma = e(i);
    %subplot(3,6,i);
    x = ijspacb(:,2+scol); y = ijspacm(:,2+scol); a = 0; b = 1.6;
    %y(x>b*std(x)) = a; x(x>b*std(x)) = a; x(y>b*std(y)) = a; y(y>b*std(y)) = a; x = x(x~=0); y = y(y~=0); %go by std of x
    plot(x,y,'.'); axis('tight');a = xlim; b = ylim; xlim([min(a(1),b(1)), max(a(2),b(2))]); ylim(xlim)   %'mo','markersize',16,'linewidth',lfs
    hold on;
    title(sprintf('sigma %.2f',round(sigma,2))); xlabel('bef'); ylabel('during');
    p = polyfit(x,y,1);
    yfit = polyval(p,x);yfit =  p(1) * x + p(2);
    plot(x,yfit)
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq = 1 - SSresid/SStotal;    
    %legend((sprintf('r = %.1f %.1f', polyfit( x, y,1) )  ),'Location','best');
    legend( sprintf('r = %.1f res = %.2f',p(1),rsq)  ,'Location','best');
    hold off;%axis('equal');%axis('square');
    
end
suptitle('spatial corr at 0,0 bef vs dur by gauss smoothing sigma');

 %plot(size(cc,2)/2,size(cc,1)/2,'md','MarkerFaceColor','w','MarkerSize',7);
%imagesc(az, imgaussfilt(c1.rm,1)); %RM


% ijcb = []; ijcm = [];
%     for ri = 1:283
%         for ci = 1:287
%             if ~isempty(corrBefore{ri,ci})
%                 t = corrBefore{ri,ci};
%                 ijcb(end+1,:) = [ri,ci,t];
%                 t = corrMid{ri,ci}; 
%                 ijcm(end+1,:) = [ri,ci,t];
%                 %ijcbmpbm(end+1,:) = [ri ci cb cm pb pm ];
%             end
%         end
%     end