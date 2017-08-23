function   shufflingMain()
    warning off;
    %%%
    binsize = 45;
    %fn = sprintf('C:\\Noam\\Data\\muscimol\\noam\\cells_%dmin_c_midall_gridscore.mat',binsize); first set
    fn = sprintf('C:\\Noam\\Data\\muscimol\\noam\\cells_%dmin_d_patchtraj_rayleigh',45);
    load(fn);
    [groups ~] = findSimultaneouslyRecordedCells(cells);
    gridThreshBef = 0.3; gridThreshMid = 0.25;
    
    %[gids, cids] = groupsByGridThresh(groups, gridThreshBef, gridThreshMid); % >grid, <mid
    
    iCbmIuRbm250u = {}; 
    for i = 1:length(groups)
        i
        g = groups{i}; %AM I ITERATIONG THROUGH ALL?
        for j = 1:length(g)-1
            for k = j+1:length(g)
                cb = shuffleTimeCorrelations (g(j).before,g(k).before, 1, 200, 0.006);
                cm = shuffleTimeCorrelations (g(j).midall,g(k).midall, 1, 200, 0.006);
                iu = Intersection_Over_Union(g(j).before.module.x, g(j).before.module.y, g(k).before.module.x,g(k).before.module.y);
                iCbmIuRbm250u{g(j).ind,g(k).ind} = [g(j).ind g(k).ind cb cm iu... 
                    g(j).before.rayleigh_score g(k).before.rayleigh_score g(j).midall.rayleigh_score g(k).midall.rayleigh_score];
            end
        end
    end
    
     iud = []; iun = []; iuo = []; ium = []; iua = [];                                       
     for i = 1:length(groups)
        g = groups{i}; %AM I ITERATIONG THROUGH ALL?
        for j = 1:length(g)-1
            for k = j+1:length(g)
                i1 = g(j).ind;  i2 = g(k).ind; gm1 =  g(j).midall.gridscore; gm2 =  g(k).midall.gridscore;
                c1 = g(j).before; c2 = g(k).before; cb = iCbmIuRbm250u{i1,i2}(3); cm = iCbmIuRbm250u{i1,i2}(4);
                iu = Intersection_Over_Union(c1.module.x,c1.module.y, c2.module.x,c2.module.y);
                %t = [iu i1 i2 cb cm];
                t = iCbmIuRbm250u{i1,i2};
                if      c1.gridscore >= gridThreshBef && gm1 <= gridThreshMid && c2.gridscore >= gridThreshBef  && gm2 <= gridThreshMid
                    iud(end +1,:) = t;
                elseif  c1.gridscore >= gridThreshBef && gm1 >  gridThreshMid && c2.gridscore >= gridThreshBef  && gm2 >  gridThreshMid
                    iun(end +1,:) = t;
                elseif (c1.gridscore >= gridThreshBef && gm1 <= gridThreshMid && c2.gridscore >= gridThreshBef  && gm2 >  gridThreshMid) ||...
                       (c1.gridscore >= gridThreshBef && gm1 >  gridThreshMid && c2.gridscore >= gridThreshBef  && gm2 <  gridThreshMid)
                    iuo(end +1,:) = t;
                else
                    ium(end +1,:) = t;
                end
                iua(end +1,:) = t;
            end
        end
     end
    %iCbmIuRbm250u              i1 i2 cb cm iu r1b r2b r1m r2m 
    %                           1  2  3  4  5  6   7   8   9
    i1i = 1; i2i = 2; cbi = 3; cmi = 4; iui = 5; r1bi = 6; r2bi = 7; r1mi = 8; r2mi = 9;
    
    % RAYLEIGH + MODULE
    criuBef = {}; criuMid = {};
    criuB =[]; criuM = [];
    for i = 1:length(gids)
        gis = gids(i)
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
    row = 4; col = 3; ii =1; figure; a = [];
    
    t = iud; lims = [min(min(t(:,cbi), min(t(:,cmi)))) max(max(t(:,cbi)), max(t(:,cmi)))] ;
    tg = t(t(:,r2mi)>=0.5 & t(:,r1mi)>=0.5,:); tl = t(t(:,r2mi)<0.5 | t(:,r1mi)<0.5,:);
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(t(:,r1mi),t(:,r2mi),'ro'); 
    xlim([0 1]); ylim([0 1]); title(sprintf('decreasing gs pairs (%.2f>%.2f), rayleigh i1 vs i2 MID',gridThreshBef,gridThreshMid));  
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(tg(:,cbi),tg(:,cmi),'ro'); title('rayleigh > 0.5 corr b vs m');
    hold on; plot([min(tg(:,cbi)) max(tg(:,cbi))], polyval( polyfit( tg(:,cbi),tg(:,cmi),1),[min(tg(:,cbi)) max(tg(:,cbi))] ) );
    xlim(lims); ylim(lims);
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(tl(:,cbi),tl(:,cmi),'ro'); title('rayleigh < 0.5 corr b vs m');
    hold on; plot([min(tl(:,cbi)) max(tl(:,cbi))], polyval( polyfit( tl(:,cbi),tl(:,cmi),1),[min(tl(:,cbi)) max(tl(:,cbi))] ) ); 
    xlim(lims); ylim(lims);
    
    t = iua; lims = [min(min(t(:,cbi), min(t(:,cmi)))) max(max(t(:,cbi)), max(t(:,cmi)))] ;
    tg = t(t(:,r2mi)>=0.5 & t(:,r1mi)>=0.5,:); tl = t(t(:,r2mi)<0.5 | t(:,r1mi)<0.5,:);
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(t(:,r1mi),t(:,r2mi),'ro'); title('all pairs rayleigh i1 vs i2 MID'); 
    xlim([0 1]); ylim([0 1]); 
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(tg(:,cbi),tg(:,cmi),'ro'); title('rayleigh > 0.5 corr b vs m');
    hold on; plot([min(tg(:,cbi)) max(tg(:,cbi))], polyval( polyfit( tg(:,cbi),tg(:,cmi),1),[min(tg(:,cbi)) max(tg(:,cbi))] ) );
    xlim(lims); ylim(lims);
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(tl(:,cbi),tl(:,cmi),'ro'); title('rayleigh < 0.5 corr b vs m');
    hold on; plot([min(tl(:,cbi)) max(tl(:,cbi))], polyval( polyfit( tl(:,cbi),tl(:,cmi),1),[min(tl(:,cbi)) max(tl(:,cbi))] ) ); 
    xlim(lims); ylim(lims);
    
        t = iud; lims = [min(min(t(:,cbi), min(t(:,cmi)))) max(max(t(:,cbi)), max(t(:,cmi)))] ;
    tg = t(t(:,r2mi)>=0.5 & t(:,r1mi)>=0.5,:); tl = t(t(:,r2mi)<0.5 | t(:,r1mi)<0.5,:);
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(t(:,r1bi),t(:,r2bi),'bo'); title('decreasing pairs, rayleigh i1 vs i2 BEF'); 
    xlim([0 1]); ylim([0 1]); 
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(tg(:,cbi),tg(:,cmi),'bo'); title('rayleigh > 0.5 corr b vs m');
    hold on; plot([min(tg(:,cbi)) max(tg(:,cbi))], polyval( polyfit( tg(:,cbi),tg(:,cmi),1),[min(tg(:,cbi)) max(tg(:,cbi))] ) );
    xlim(lims); ylim(lims);
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(tl(:,cbi),tl(:,cmi),'bo'); title('rayleigh < 0.5 corr b vs m');
    hold on; plot([min(tl(:,cbi)) max(tl(:,cbi))], polyval( polyfit( tl(:,cbi),tl(:,cmi),1),[min(tl(:,cbi)) max(tl(:,cbi))] ) ); 
    xlim(lims); ylim(lims);
    
    t = iua; lims = [min(min(t(:,cbi), min(t(:,cmi)))) max(max(t(:,cbi)), max(t(:,cmi)))] ;
    tg = t(t(:,r2mi)>=0.5 & t(:,r1mi)>=0.5,:); tl = t(t(:,r2mi)<0.5 | t(:,r1mi)<0.5,:);
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(t(:,r1bi),t(:,r2bi),'bo'); title('all pairs rayleigh i1 vs i2 BEF'); 
    xlim([0 1]); ylim([0 1]); 
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(tg(:,cbi),tg(:,cmi),'bo'); title('rayleigh > 0.5 corr b vs m');
    hold on; plot([min(tg(:,cbi)) max(tg(:,cbi))], polyval( polyfit( tg(:,cbi),tg(:,cmi),1),[min(tg(:,cbi)) max(tg(:,cbi))] ) );
    xlim(lims); ylim(lims);
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(tl(:,cbi),tl(:,cmi),'bo'); title('rayleigh < 0.5 corr b vs m');
    hold on; plot([min(tl(:,cbi)) max(tl(:,cbi))], polyval( polyfit( tl(:,cbi),tl(:,cmi),1),[min(tl(:,cbi)) max(tl(:,cbi))] ) ); 
    xlim(lims); ylim(lims);
    
    
    
    %INTERSECTION UNION
    row = 5; col = 3; ii = 1; 
    figure; 
     t = iud; t = t(t(:,iui)>0,:); tg = t(t(:,iui)>=0.6,:); tl = t(t(:,iui)<0.6,:); 
    a(end+1) = subplot(row,col,ii);ii=ii+1; hist(t(:,iui),20); 
    title(sprintf('decreasing gs pairs (%.2f>%.2f) i/u',gridThreshBef,gridThreshMid)); xlim([0 1]);
    a(end+1) = subplot(row,col,ii);ii=ii+1;plot(tg(:,cbi),tg(:,cmi),'ro');title('iu >= 0.6 corr b vs m');
    hold on; plot([min(tg(:,cbi)) max(tg(:,cbi))], polyval( polyfit( tg(:,cbi),tg(:,cmi),1),[min(tg(:,cbi)) max(tg(:,cbi))] ) );
    axis equal;
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(tl(:,cbi),tl(:,cmi),'bo'); title('iu < 0.6 corr b vs m');
    hold on; plot([min(tl(:,cbi)) max(tl(:,cbi))], polyval( polyfit( tl(:,cbi),tl(:,cmi),1),[min(tl(:,cbi)) max(tl(:,cbi))] ) ); 
    axis equal;
     t = iun;t = t(t(:,iui)>0,:); tg = t(t(:,iui)>=0.6,:); tl = t(t(:,iui)<0.6,:); 
    a(end+1) = subplot(row,col,ii);ii=ii+1; hist(t(:,iui),20); title('non decreasing pairs i/u');xlim([0 1]);
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(tg(:,cbi),tg(:,cmi),'ro'); title('iu >= 0.6 corr b vs m');
    hold on; plot([min(tg(:,cbi)) max(tg(:,cbi))], polyval( polyfit( tg(:,cbi),tg(:,cmi),1),[min(tg(:,cbi)) max(tg(:,cbi))] ) );
    axis equal;
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(tl(:,cbi),tl(:,cmi),'bo'); title('iu < 0.6 corr b vs m');
    hold on; plot([min(tl(:,cbi)) max(tl(:,cbi))], polyval( polyfit( tl(:,cbi),tl(:,cmi),1),[min(tl(:,cbi)) max(tl(:,cbi))] ) );  
    axis equal;
     t = iuo; t = t(t(:,iui)>0,:); tg = t(t(:,iui)>=0.6,:); tl = t(t(:,iui)<0.6,:); 
    a(end+1) = subplot(row,col,ii);ii=ii+1; hist(t(:,iui),20); title('one decreasing pairs i/u');xlim([0 1]);
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(tg(:,cbi),tg(:,cmi),'ro'); title('iu >= 0.6 corr b vs m');
    hold on; plot([min(tg(:,cbi)) max(tg(:,cbi))], polyval( polyfit( tg(:,cbi),tg(:,cmi),1),[min(tg(:,cbi)) max(tg(:,cbi))] ) );
    axis equal;
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(tl(:,cbi),tl(:,cmi),'bo'); title('iu < 0.6 corr b vs m');
    hold on; plot([min(tl(:,cbi)) max(tl(:,cbi))], polyval( polyfit( tl(:,cbi),tl(:,cmi),1),[min(tl(:,cbi)) max(tl(:,cbi))] ) );
    axis equal;
     t = ium; t = t(t(:,iui)>0,:); tg = t(t(:,iui)>=0.6,:); tl = t(t(:,iui)<0.6,:); 
    a(end+1) = subplot(row,col,ii);ii=ii+1; hist(t(:,iui),20); title('rest of pairs i/u');xlim([0 1]);
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(tg(:,cbi),tg(:,cmi),'ro'); title('iu >= 0.6 corr b vs m');
    hold on; plot([min(tg(:,cbi)) max(tg(:,cbi))], polyval( polyfit( tg(:,cbi),tg(:,cmi),1),[min(tg(:,cbi)) max(tg(:,cbi))] ) );   
    axis equal;
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(tl(:,cbi),tl(:,cmi),'bo'); title('iu < 0.6 corr b vs m');
    hold on; plot([min(tl(:,cbi)) max(tl(:,cbi))], polyval( polyfit( tl(:,cbi),tl(:,cmi),1),[min(tl(:,cbi)) max(tl(:,cbi))] ) );
    axis equal;
     t = iua; t = t(t(:,iui)>0,:); tg = t(t(:,iui)>=0.6,:); tl = t(t(:,iui)<0.6,:); 
    a(end+1) = subplot(row,col,ii);ii=ii+1; hist(t(:,iui),20); title('all pairs i/u');xlim([0 1]);
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(tg(:,cbi),tg(:,cmi),'ro'); title('iu >= 0.6 corr b vs m');
    hold on; plot([min(tg(:,cbi)) max(tg(:,cbi))], polyval( polyfit( tg(:,cbi),tg(:,cmi),1),[min(tg(:,cbi)) max(tg(:,cbi))] ) );
    axis equal;
    a(end+1) = subplot(row,col,ii);         plot(tl(:,cbi),tl(:,cmi),'bo'); title('iu < 0.6 corr b vs m');
    hold on; plot([min(tl(:,cbi)) max(tl(:,cbi))], polyval( polyfit( tl(:,cbi),tl(:,cmi),1),[min(tl(:,cbi)) max(tl(:,cbi))] ) );
    axis equal;
    
    
    
    t = findObj(gcf,'type','axes');
    for i = 1:length(t)
        a = t(i);
        if mod(i-1,3) == 0
            mod(i-1,3)
            axis(a, 'equal');
        end
    end
    
    
    
    figure; %hist iu
    subplot(1,3,1);hist(iud,20); title('IU pairs decreasing'); xlabel('intersection / union');
    subplot(1,3,2);hist(iun(iun>0),20); title('IU pairs non-decreasing');
    subplot(1,3,3);hist(iua(iua>0),20); title('IU pairs all'); 
    
    % SHUFFLING
    gcorrBefore = {}; gcorrMid = {};
     corrBefore = {};  corrMid = {};
    n = 250;
    tic
    for i = 1:length(gids)
        gis = gids(i)
        cis = cids{gis}
        g = groups{gis};
        g = g(cis);
        for j = 1:length(g)-1
            for k = j+1:length(g)
                [i j k] 
                gcorrBefore{k,j,i} = shuffleTimeCorrelations (g(j).before, g(k).before, n, 200, 0.006);
                gcorrMid{k,j,i}    = shuffleTimeCorrelations (g(j).midall, g(k).midall, n, 200, 0.006);
                corrBefore{g(j).ind,g(k).ind} = gcorrBefore{k,j,i};
                   corrMid{g(j).ind,g(k).ind} = gcorrMid{k,j,i};
            end
        end
    end
    disp('************************');
    toc
      
    
%     t = corrBefore{11,12};
%     [s,I] = sort(abs(t),'descend');
%     p = round(find(I==1)/length(I),2) %index of non shuffled 
    i = 11, j = 12, c1 = g(1).before, c2 = g(2).before, lag = 200, sigma = 0.006
    [r, c] = size(corrBefore);
    pvalsBefMid = cell(r,c); pbs = []; pms = [];
    for i = 1:r
        for j = 1:c
            if ~isempty(corrBefore{i,j})
                t = corrBefore{i,j};
                [s,I] = sort(abs(t),'descend');
                pb = round(find(I==1)/length(I),2); %index of non shuffled
                t = corrMid{i,j};
                [s,I] = sort(abs(t),'descend');
                pm = round(find(I==1)/length(I),2);
                pvalsBefMid{i,j} = [pb pm]; pbs(end+1) = pb; pms(end+1) = pm;
            end
        end
    end
    figure;
    subplot(2,1,1);hist(pbs*100,100); title('p shuffle before 250 update');
    subplot(2,1,2);hist(pms*100,100); title('p shuffle mid 250 update');
    
    disp('');
    %}
end

function pcor0 = shuffleTimeCorrelations (c1, c2, n, lag, sigma)
    %n=100;
    %lags = 0.2; sigma = 0.006; %lag in s; determined empirically
    %convert to ms and start at 1
    %c2 = cells{71}.before; c1 = cells{75}.before; %BEFORE   
    s1=floor(c1.st*1000)+1; s2=floor(c2.st*1000)+1; 
    mint=min(min(s1),min(s2));
    s1 = s1-mint+10; s2 = s2-mint+10; %pad edges
    maxt=max(max(s1),max(s2))+10; %MAKE CORRELATIONS SMALLER? pad ends
    %remove overlapping spikes
    [s1i s2i] = removeOverlappingSpikes(s1,s2, 1); s1=s1(s1i); s2=s2(s2i);
    %create train for c1
    train1 = zeros(1,maxt); train1(s1)=1;
    train2 = zeros(1,maxt); train2(s2)=1;
    %delta to shift train by
    d = floor(length(train1)/n) ;%right?  %amount to shift by each step
    %trains=zeros(n+1,maxt); trains(i,:)=train2shifted; save all
    pcor0 = []; 
    %t2shifted = zeros(n+1,maxt);
    %win=hamming(15);
    for i = 1:n+1 %n+1 because wi = (i-1)*d; to include 0 shift
        %t2shifted(i,:) = conv(train2([(i-1)*d+1:length(train1) 1:(i-1)*d]),win,'same');
    end 
    %pcor=xcov(train1Ham, train2Ham,lag,'coef');
    
    for i = 1:n %n+1 because wi = (i-1)*d; to include 0 shift
        sprintf('done %.2f\%',round(i/(n+1),2));
        wi = (i-1)*d; %first one calculated is 0
        shiftedi = [wi+1:length(train1) 1:wi]; %cyclic
        train2shifted = train2(shiftedi); %trains(i,:)=train2shifted;
        pcor = timeCorrelationSmoothed(train1,train2shifted,15,lag,sigma);
        pcor0(i) = pcor(round(length(pcor)/2));
        %{
        %smoothing spike train
        win=hamming(15);
        %train1Ham = train1; train2ham = train2;
        train1Ham=conv(train1,win,'same');train2Ham=conv(train2shifted,win,'same'); %win = win/sum(win);
        %  figure, plot(train1Ham), hold on, plot(train1 + 1);
        %correlation
        %pcor=xcov(train1Ham, train2Ham,lag,'coef'); % TRY AS MATRIX
        pcor=xcorr(train1Ham, train2Ham,lag,'coef');
        %p{i,3} = [pmx (pmxi - ((length(p{i,1})-1)/2) - 1)/1000]; % dont ask
        %smoothing correlation
        [b,a] = butter(6,sigma); %LOW PASS 0.15%6th order, fc/fs/2 determined empiracally 
        pcor = filtfilt(b,a,pcor);
        
        pcor0(i) = pcor(round(length(pcor)/2)); %corr at 0
        %count_xcor{i,1}=xcorr(train1,train2, round(1000*lags))'; %how many times spike at same point
        %}
    end
    %disp(pcor0);
    %[mx, imx] = max(abs(pcor0))
    %time_scale=( (1:length(p{i,1})) - ((length(p{i,1})-1)/2) - 1)/1000; %goes from -lag 0 +lag in ms
    %figure;
     %plotmatrix(time_scale',cell2mat(p(:,1))','-');
     %figure
     %plotmatrix((1:maxt)',trains','-');
    %plotmatrix(cell2mat(p(1,2))',cell2mat(p(:,1))')
end


function [gids, cids] = groupsByGridThresh(groups, gridThreshBef, gridThreshMid)
    gids = []; cids = {};
    for i = 1:length(groups)
        g = groups{i};
        t = [];
        for j = 1:length(g)
            if g(j).before.gridscore > gridThreshBef...
            && g(j).midall.gridscore < gridThreshMid
            %if g(j).before.gridscore < g(j).midall.gridscore 
                t(end+1) = j;
            end
        end
        if length(t) >= 2 %only groups with at least 2
            cids{i} = t; 
            gids(end+1) = i;
        end
    end
end