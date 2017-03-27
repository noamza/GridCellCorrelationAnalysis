%%%%%%%%%
% N Z A %
% grid cells: 228
% pairs: 749
%%% DATASTATS() !!
%%%%%%%%%

%{
** TO DO

- fix extrema2

    ****Questions***
- Nan errors in gridscore, etc
- Method for correlations that don't fit

* Check if firing rate shuts down by location of environment in grid cell

% differences: mid length
               after existence
               after lengths
               area
               number of middle fields
               g_12048_030708 * 54956 31508 54956 * after
% same
               no overlapping tet/ch combo (up to 9?)
               before length same
               
read/study:
        correlation + auto + cross, covatiance.
        Bonavie paper   < < < < <   
%}

function musciomol()
    warning off;
    params.dir_load =...
        'C:\Noam\Data\muscimol\DB_musc_MEC\';
    params.dir_save =....
        'C:\Noam\Data\muscimol\noam_output\';
    
    binsize = 10;                            % B I N S I Z E
    
    fn = sprintf('C:\\Noam\\Data\\muscimol\\noam\\cells_%dmin_a.mat',binsize);
    %    tic %1800 sec   %LOAD FROM FILES
    %{
    tic;
    files = dir(strcat(params.dir_load,'DB*.mat'));
    for i= 1:length(files)
        data{i} = load(strcat(params.dir_load,files(i).name));
        data{i}.db.ind = i;
        cells{i} = process(data{i}.db, binsize); 
        fprintf('%d %.1f \n',i, (i*100/length(files)));
    end
    save(fn,'cells');
    toc
    %}
    
    %fprintf('loading **bin size**%d*\n',binsize); %ascii 48
    %tic; cells = load(fn); cells = cells.cells; toc;
    
    %saveSpikesByDirection(params, 8);
    
    plotByDirectionMain();
    

    stop
    groups = find_simultaneously_recorded_cells(cells);
    k = fieldnames(groups);
   
    disp('done and done'); 
    
    %TEMPORAL CORRELATION
    %{
    figure();
    tic
    for i = 1:length(k)
        %disp(i);
        colormap jet;
        g = groups.(k{i});
        [mc, lc] = plotGroupTemporralCorrs(g, 500);
%         subplot(10,9,2*i-1);
%         imagesc(lc./max(lc(:))); axis off; axis equal; caxis([-1 1]);
%         title(sprintf('i%d:[%.1f %.1f]',i, min(lc(:)), max(lc(:))));
%         set(gca,'FontSize',8)
%         subplot(10,9,2*i);
%         imagesc(mc); axis off; axis equal; %caxis([-1 1]);
%         [m, ind] = max(mc(:));
%         title(sprintf('i%d:[%.1f %.3f]',i, lc(ind), m));
%         set(gca,'FontSize',8)
         subplot(5,8,i);
         imagesc(mc,[0,0.5]); axis off; axis equal; %caxis([-1 1]);
         [m, ind] = max(mc(:));
         title(sprintf('i%d(%.3f)(%.1f)',i, m, lc(ind)));
         set(gca,'FontSize',9); 
        
        fprintf('%d %.1f %.3f\n',i,lc(ind),m);
    end
    set(suptitle('Highest Temporal-Correlation/Lag for Each Pair in Group, Normalized'),'Interpreter', 'none');
    toc; disp('time corr');
    %}
    
    %find width of plot
    midmax = 0;
    for i = 1:length(k)
        g = groups.(k{i});
        g = sortByMiddle(g);
        length(g(1).middle);
        midmax = max(midmax ,length(g(1).middle)); 
    end
    %disp(midmax); %lots of 8
    %%%%
    midmax = 6;%%%%%%%
    
    %i = 5; plot_group(groups.(k{i}), k{i}, i, binsize, midmax);
    
    %plot groups
    for i = 1:length(k)
        fprintf('%s %d %d\n', k{i}, i, length(groups.(k{i})));
        plot_group(groups.(k{i}), k{i}, i, binsize, midmax);
    end

    %print group properties
    for i = 5%1:length(k);
        g = groups.(k{i});
        t = g(1); t = length(t.middle); %t = length(g);
        %t.a; length(t.after.px) .exists length(t.middle)
        fprintf('%s * %d \n',k{i} ,t);
        %fprintf('%s %d', k{i}, length(g));
        for j = 1:length(g) %in group
            tt = g(j); tt = length(tt.middle);
            %.a length(tt.after.px) .exists length(tt.middle)
            if t ~= tt               %not(strcmp(t,tt))
                fprintf('%d ', tt);
                t = tt;
            end; %fprintf('%d %d\n', g(j).tet, g(j).cel);
        end; disp('*'); 
    end;



    %pairs
    tic
    pairs = find_pairs(cells);%find_pairs_of_cells_mosimol(data);
    d = 47; %65 23 86 235
    for i = 48:49 %47:length(pairs)
        r1i = pairs(i,1); r2i = pairs(i,2);
        fprintf('%.1f of %d ',i*100/length(pairs),length(pairs));
        r1 = cells{r1i};%process(data{r1i}.db);
        r2 = cells{r2i};%process(data{r2i}.db);
        %plot_pair(r1, r2, i);
    end
    toc
    disp('muscimol done');
    %{
        - Gilad's paper
        - to what extent do cells disappear together
        - max simul recorded
        - module
    %}
end

function saveSpikesByDirection(params, degbins)
    tic; files = dir(strcat(params.dir_load,'DB*.mat'));
    prnt = true;
    for i= 1:length(files)
        fprintf('%d %.1f \n',i, (i*100/length(files)));
        a = load(strcat(params.dir_load,files(i).name));
        %a = load(strcat(params.dir_load,'DB_MUSC_MEC_001_11468_010306_t2c1.mat '));
        p = a.db.B(1).pos_data; s = a.db.B(1).spike_data;

        %HEAD DIRECTION &&&&&&
        % (pos_t,pos_x,pos_y,pos_x2,pos_y2,spk_x,spk_y,spk_x2,spk_y2, bins)
        [HD_rate_bins,HD_count, HD_time, ang_bins]= computeHeadDirectionality...
            (p.t,p.x1,p.y1,p.x2,p.y2,s.x,s.y,s.x2,s.y2, degbins);

        %MOVING DIRECTION
        assert(length(s.x) == length(s.ts), 'length(s.x) == length(s.ts)');
        a = rat_trajectory(p.x1, p.x2, p.y1, p.y2, p.t, s.x, s.x2, s.y, s.y2, s.ts);
        %tic; cells = load(fn); cells = cells.cells; toc; a = cells{1}.before;
        px = double(a.px); py = double(a.py); pt = a.pt; st = a.st;
        % computing velocity
        parms.smoothing_parameter_for_velocity = 10e-3; %10e-5;
        tx = csaps( 1:length(pt),px);%parms.smoothing_parameter_for_velocity );
        ty = csaps( 1:length(pt),py);%parms.smoothing_parameter_for_velocity ); %param = 10e-5 or 10e-3
        %spline and smoothing, off by one binning and smoothing
        % Compute the velocity via a derivative of the pp-form of the spline (fnder.m):
        dt=median(diff(pt));
        vx = fnval( fnder(tx),1:length(pt))/dt; %tx
        vy = fnval( fnder(ty),1:length(pt))/dt; %ty
        %vx = diff(px)/dt; vx = [0 vx];
        %vy = diff(py)/dt; vy = [0 vy];
        bins = 8; parms.num_of_direction_bins = bins;
        [MD_rate_bins,MD_count,MD_time, ang_bins]=...
            computeMovingDirectionality(pt,px,py,st,vx,vy,parms); %Compute_Moving_Directionality
        
        %correct
        phd = atan2(p.y2-p.y1,p.x2-p.x1);
        pmd = atan2(vy,vx)';
        cellByDeg{i}.ind = i; cellByDeg{i}.ang_bins = ang_bins;
        cellByDeg{i}.HD_rate = HD_rate_bins;  cellByDeg{i}.MD_rate = MD_rate_bins;
        if prnt
            h = figure('Position', [0, 0, 1000, 500]);%(2*m+5)*120 , 140*r + height]); %
            set(gca,'LooseInset', get(gca,'TightInset'));colormap jet;
            titl = sprintf('Head_Moving_Direction_%dbins_i%d', bins, i);
            subplot(2,3,1);bar(ang_bins,HD_rate_bins); title('HD firing rate');
            subplot(2,3,2);bar(ang_bins,HD_time);       title('Time');
            subplot(2,3,3);bar(ang_bins,HD_count);      title('Spikes');     
            subplot(2,3,4);bar(ang_bins,MD_rate_bins); title('MD firing rate');
            subplot(2,3,5);bar(ang_bins,MD_time);       title('Time');
            subplot(2,3,6);bar(ang_bins,MD_count);      title('Spikes');        
            ax = findobj(gcf,'Type','Axes');
            for q=1:length(ax)
                set(ax(q),'FontSize',9);%axis(ax(q),'equal');%axis(ax(q),'off');%title(ax(i),{'Very Nice'})
            end
            set(suptitle(titl),'Interpreter', 'none');
            sdir = 'C:\Noam\Output\muscimol\HDMD\'; debug = '';
            filename = sprintf('%s%s%s.png',sdir, debug, titl); disp(filename);
            h.PaperPositionMode = 'auto'; print(filename, '-dpng','-r0');
            close(h)
        end
    end
    %SAVE VAR
    fn = sprintf('C:\\Noam\\Data\\muscimol\\noam\\direction_spikes_%ddegs.mat',degbins);
    disp(fn);
    save(fn,'cellByDeg');
end
function correlation_sanity_check(c1, c2)
    nbins = 4;
    %MOVING DIRECTION
    smooth = 10e-7;%10e-3 %10e-5; %Empiracally so bins don't change rapidly for 8 bins!
    tx = csaps( 1:length(c1.pt),double(c1.px), smooth);%10e-3; %10e-5;
    ty = csaps( 1:length(c1.pt),double(c1.py), smooth);
    dt = median(diff(c1.pt));
    vx = fnval( fnder(tx),1:length(c1.pt))/dt; 
    vy = fnval( fnder(ty),1:length(c1.pt))/dt;
    pos_md = atan2d(vy,vx);%wrapTo360(rad2deg(atan2d(vy,vx)));    
    %TO 360 
    pos_md = wrapTo360(rad2deg(pos_md)); %check??
    %BIN BY DIRECTIONS
    bin_edges = (0:360/nbins:360);% - 180/nbins; bin_edges(1) = 0; bin_edges(end +1) = 360; 
    pos_by_bin = discretize(pos_md, bin_edges)';% pos_by_bin(pos_by_bin==nbins+1) = 1; %conbines two last bins into one, to wrap pi
    a = diff(c1.pt);
    unique(a)
    
    
    
    edges = round(edges);
    direction_array = interp1(c1.pt, unwrap(pos_md'),edges);
    
    %SPIKES
    c1.st = c1.st(c1.st>=min(c1.pt)&c1.st<=max(c1.pt)); %removes spikes out of bound pos times
    c2.st = c2.st(c2.st>=min(c1.pt)&c2.st<=max(c1.pt)); %removes spikes out of bound pos times
    %to 360
    spk1_md = interp1(c1.pt, unwrap(pos_md'),c1.st);
    spk2_md = interp1(c1.pt, unwrap(pos_md'),c2.st);
    spk1_md = wrapTo360(rad2deg(spk1_md)); spk2_md = wrapTo360(rad2deg(spk2_md));
    %FIRING RATE
    spkt1 = c1.st; spkt2 = c2.st;
    spkt1 = floor(spkt1*1000)+1; spkt2=floor(spkt2*1000)+1; %in MILLISECS
    max_time=max(max(spkt1),max(spkt2));
    %firing rate;
    bin_msecs = 100;
    firing1 = zeros(max_time,1); firing2 = zeros(max_time,1);
    firing1(spkt1) = spkt1;  firing2(spkt2) = spkt2;    
    [firing1, edges] = histcounts(firing1,round(max_time/bin_msecs));%(x,nbins) 
    firing2 = histcounts(firing2,round(max_time/bin_msecs));%now in units of binsecs
    firing1(1) = 0; firing2(1) = 0; %spikes in first bin_msecs ignored
    
    %   
    spk1_by_bin = discretize(spk1_md, bin_edges);  %spk1_by_bin(spk1_by_bin==nbins+1) = 1;
    spk2_by_bin = discretize(spk2_md, bin_edges);    
    
    [P,F] = pwelch(y);
    helperFilterIntroductionPlot1(F,P,[60 60],[-9.365 -9.365],...
  {'Original signal power spectrum', '60 Hz Tone'})
    fc = 500;
    fs = 10000;
    figure;
    plot(y); hold on;
    plot(filter(b,a,t));
    [z,p,k] = butter(10,0.05);
    sos = zp2sos(z,p,k);
    grpdelay(sos,128);
    
    
end


%{
    % correlation bt 2 simultaneously recorded cells, same module(?) not overlapping(?)
    % for each 8 directions, temporal correlation between the 2 cells, vary amount of time
    % use moving direction
    % correction not for this
    % correction cell by cell
 - Take group
 - Find grid cells
 - Take one cell, for each direction, try to correlate each other cell
 - for each group x for each cell x for each direction x all other cells
 - for g x c x d x cs
 - COMPARE: 75 76, 226 228 
%}
function plotByDirectionMain()    
    fn = sprintf('C:\\Noam\\Data\\muscimol\\noam\\cells_%dmin_a.mat', 10);
    fprintf('loading **bin size**%d*\n',10); %ascii 48
    tic; cells = load(fn); cells = cells.cells; toc;
%     correlation_sanity_check(cells{75}.before, cells{76}.before);
    groups = find_simultaneously_recorded_cells(cells);
    k = fieldnames(groups);
    params.time_bin_secs = 0.01;
    params.lag_max_secs = 10;
    nbins = 4; params.number_degree_bins = nbins;
    %load cells
    %%%% for each group x for each cell x for each direction x all other cells
    for l = 1:length(k)
        g = groups.(k{l});
        %fprintf('%s %d %d\n', k{l}, l, length(groups.(k{l})));
        prnt = true;
        rmt = 1.5; % Rate Map Threshold Hz
        gridThresh =  0.5; % was .7
        good = [g(1)]; bad = [g(1)]; %find low gridscore cells
        for j = 1:length(g) %put low grid scores at end
            if g(j).before.gridscore > gridThresh
                good(end+1) = g(j);
            else bad(end+1) = g(j);end; end
        bad = bad(2:end); good = good(2:end); %[good(2:end) bad(2:end)]; Only show good cells  
        %FIGURE
        %%% plot each cell %%%
        r = length(good)-1; c = r; if r == 1; r=2; end; 
        prnt = 0;
        if ~isempty(good) && length(good) >= 2
            prnt = 1;
            h = figure('Position', [0, 0, 2000, 1000]);%(2*m+5)*120 , 140*r + height]); %
            set(gca,'LooseInset', get(gca,'TightInset'));
            colormap jet;
            %for each cell, for each direction x each other cell
            for i = 1:length(good)-1
                for j = i+1:length(good)
                    %i; j;
                    [p co] = compareByMovingDirection(good(i).before, ...
                        good(j).before, params); %j

                    subplot(r, c, c*(i-1) + j-1);
                    X = [0:length(p{1})-1]'-(length(p{1})-1)/2; X=X*params.time_bin_secs;
                    Y = reshape(cell2mat(p(:,1)),[],nbins);
                    [~,ax,AX] = plotmatrix(X,Y,'-');
                    for q=1:length(ax)
                        %axis(ax(q),'off');
                        axis(ax(q),'tight')
                        ylim(ax(q),[min(Y(:)) max(Y(:))]);
                        ylabel(ax(q),sprintf('%d°',p{q,2}(1)));
                        set(ax(q),'ytick',[]);
                        set(ax(q),'yticklabel',[]);
                        set(ax(q),'FontSize',6);%axis(ax(q),'equal');%axis(ax(q),'off');%title(ax(i),{'Very Nice'})
                    end
                    if i==1 && j==2
                        axis(AX,'tight')
                        %set(AX,'FontSize',8)
                        set(AX,'fontweight','bold');
                        ylabel(AX,sprintf('MD bin (b=%.1f°)',360/nbins));
                        xlabel(AX,sprintf('%s','Xcorrelation lag (s)'));
                    end
                    title(sprintf('c%d*c%d maxxcorr=%.3f',good(i).ind,good(j).ind,max(Y(:))));
                    
                end
            end
        end
        if prnt
            titl = sprintf('moving_direction_%dbin_xcorr_i%d_%s_ncells%d_k',nbins,l,k{l},length(good));
            set(suptitle(titl),'Interpreter', 'none'); %PUT SUP TITLE AFTER ALL SUBPLOT COMMANDS << test    
            subplot(r, c, c+1);%r*(r-1)+1);
            %ratemat
            rmt = good(1).before.rm;tick = size(rmt,1);
            xlab = [sprintf('c%d ',good(1).ind)];
            for i=2:length(good)
                rmt = [rmt ones(tick,1).*max(rmt(:)) good(i).before.rm];
                xlab =[xlab; sprintf('c%d ', good(i).ind)];
            end
            colormap jet;imagesc(rmt); pbaspect([i 1 1]);
            set(gca,'xtick',[1:i]*(tick+1) - tick/2,'xticklabel',xlab,...
                'xaxislocation','top','yticklabel',[],'fontsize',12);
            
            
            sdir = 'C:\Noam\Output\muscimol\HDMD\'; debug = '';
            filename = sprintf('%s%s%s.png',sdir, debug, titl); disp(filename);
            h.PaperPositionMode = 'auto'; print(filename, '-dpng','-r0');
            close(h);
        end
    end
    stop
end

%compares moving bins of c1 to entire c2
function [pearson_xcov, count_xcor] = compareByMovingDirection(c1, c2, params)
    nbins = params.number_degree_bins;
    smooth = 10e-7;%10e-3 %10e-5; %Empiracally so bins don't change rapidly for 8 bins!
    tx = csaps( 1:length(c1.pt),double(c1.px), smooth);%10e-3; %10e-5;
    ty = csaps( 1:length(c1.pt),double(c1.py), smooth);
    dt = median(diff(c1.pt));
    vx = fnval( fnder(tx),1:length(c1.pt))/dt; 
    vy = fnval( fnder(ty),1:length(c1.pt))/dt;
    for i = 1:length(vx)
        if sqrt((vx(i)).^2 + (vy(i)).^2) <= 3
            vx(i) = nan; vy(i) = nan;
        end
    end
    %figure();plot(diff(c1.px)/dt);hold on;plot(vx)
    pos_md = atan2d(vy,vx);%wrapTo360(rad2deg(atan2d(vy,vx)));
    %if c1.st(end)>c1.pt(end);c1.st(end)=c1.pt(end);end
    %if c2.st(end)>c1.pt(end);c2.st(end)=c1.pt(end);end
    c1.st = c1.st(c1.st>=min(c1.pt)&c1.st<=max(c1.pt));
    c2.st = c2.st(c2.st>=min(c1.pt)&c2.st<=max(c1.pt));
    spk1_md = interp1(c1.pt, unwrap(pos_md'),c1.st);
    spk2_md = interp1(c1.pt, unwrap(pos_md'),c2.st);
    %TO 360 
    pos_md = wrapTo360(rad2deg(pos_md)); %check??
    spk1_md = wrapTo360(rad2deg(spk1_md)); spk2_md = wrapTo360(rad2deg(spk2_md));
    %BIN BY DIRECTIONS
    bin_edges = (0:360/nbins:360);% - 180/nbins; bin_edges(1) = 0; bin_edges(end +1) = 360; 
    pos_by_bin = discretize(pos_md, bin_edges)';% pos_by_bin(pos_by_bin==nbins+1) = 1; %conbines two last bins into one, to wrap pi
    spk1_by_bin = discretize(spk1_md, bin_edges); 
    spk2_by_bin = discretize(spk2_md, bin_edges); 
    spk1_by_bin = spk1_by_bin(~isnan(spk1_by_bin));  %remove nans
    spk2_by_bin = spk2_by_bin(~isnan(spk2_by_bin));
    assert(sum(isnan([ ... %{pos_by_bin;}% 
        spk1_by_bin; spk2_by_bin]))==0,'nan bins compareByMovingDirection()');%+sum(isnan(spk_bins))
    assert(sum(c1.pt-c2.pt)==0,'diff times compareByMovingDirection()');
    dir_bins_spk1 = cell(nbins,1);
    dir_bins_spk2 = cell(nbins,1); 
    for i = 1:length(spk1_by_bin)
        dir_bins_spk1{spk1_by_bin(i)}(end+1) = c1.st(i);
    end
    for i = 1:length(spk2_by_bin)
        dir_bins_spk2{spk2_by_bin(i)}(end+1) = c2.st(i);
    end    
    %SPIKE TRAINS
    pearson_xcov = cell(8,2); count_xcor = cell(8,2);
    for i = 1:nbins        
        [pearson_xcov{i,1}, count_xcor{i,1}] = ...
            time_correlation(dir_bins_spk1{i},dir_bins_spk2{i},params);
        pearson_xcov{i,2} = [bin_edges(i) bin_edges(i+1)];
          count_xcor{i,2} = [bin_edges(i) bin_edges(i+1)];
        %plot(pearson_xcov); hold on%figure;plot(count_xcor)
    end
    fprintf('%.2f ',round(max(max([pearson_xcov{:,1}])), 2));
    
    max(max([count_xcor{:,1}]));
        
end


%spkt1 = g(i).before.st;
function [pxcsmooth, count_xcor] = time_correlation(spkt1,spkt2, params)
    %               generating trains
    %scales times x 1000 (ms), +1 so non-zero
    bin_secs = params.time_bin_secs;
    lag = params.lag_max_secs/bin_secs;
    spkt1=floor(spkt1*1000)+1; spkt2=floor(spkt2*1000)+1; %in MILLISECS
    min_time=min(min(spkt1),min(spkt2));
    spkt1 = spkt1-min_time+1; spkt2 = spkt2-min_time+1;%normalizes
    max_time=max(max(spkt1),max(spkt2));
    %generating the spike train
    train1=zeros(1,max_time);train2 =zeros(1,max_time); 
    train1(spkt1)=1; train2(spkt2)=1;%array at these indices(time*1000) will be 1, else 0
    %[train1,train2,count_xcor,pearson_xcor] = My_Xcor(train(1,:),train(2,:));
    assert(lag>=1,'lag=>1');
    bin_msecs = bin_secs *1000;
    %binned!
    train1 = histcounts([1 spkt1 max_time],round(max_time/bin_msecs));%(x,nbins) 
    train2 = histcounts([1 spkt2 max_time],round(max_time/bin_msecs));%now in unites of binsecs
    %SUBTRACT FIRST AND LAST ADDED SPIKE
    train1(1) = train1(1)-1;train2(1) = train2(1)-1;train1(end) = train1(end)-1;train2(end) = train2(end)-1;
    win=hamming(5); %light 3, could make 5)
    train1smooth=conv(train1,win,'same');
    train2smooth=conv(train2,win,'same');
    %generating the cross corelations normelized to pearson
    pearson_xcor=xcov(train1smooth,train2smooth,lag,'coef')'; %500 parms.max_lag NEED THIS <<<
	[b,a] = butter(6,0.03); %LOW PASS 0.15%6th order, fc/fs/2 determined empiracally 
    pxcsmooth = filtfilt(b,a,pearson_xcor);
    %close all; plot(pearson_xcor); hold on; plot(pxcsmooth);plot(xcov(train1,train2,lag,'coef')'); legend('smooth','none');
    count_xcor=xcorr(train1smooth,train2smooth, lag)'; %how many times spike at same point <<<
end


function [maxcorr, lagcorr] = plotGroupTemporralCorrs(g,lag)
    maxcorr = zeros(length(g));
    lagcorr = zeros(length(g));
    for i=1:length(g) 
        t1 = g(i).before.st;
        for j=i+1:length(g)
            t2 = g(j).before.st;
            [p, c] = time_correlation(t1,t2);
            %A(i,j) = corr(t1, t2);           
            [m, ind] = max(p);
            maxcorr(i,j) = m;
            lagcorr(i,j) = ind - (lag + 1); %+1 so that 1 - (500 + 1) = -500
        end
    end  
end


function [fig]= plotAcorrModule(c, fig, sub, s)%m,n,l
    subplot(sub(1), sub(2), sub(3)); %(m,n,l);
    imagesc(c.ac');  hold on; xlim(xlim); ylim(ylim);%axis('xy'); axis ij; axis equal; axis off;
    x=xlim; y=ylim;
    title(sprintf('%s%.2f(%d)',s, c.gridscore, length(c.sx)));
    colormap jet;
    lw = 0.5;
    major = c.module.major_ax; minor = c.module.minor_ax; hex_peaks = c.module.hex_peaks;
    phi = c.module.angle; %+ pi/2;
    co = 'k';
    %center point
    if hex_peaks ~= -1 %CHANGE TO EXISTS
        x0=hex_peaks(7,1); y0=hex_peaks(7,2);
        beta = phi; sinbeta = sin(beta); cosbeta = cos(beta);
        alpha =0: pi/100:2*pi; sinalpha = sin(alpha);cosalpha = cos(alpha);
        x1 = x0 + (major * cosalpha * cosbeta - minor * sinalpha * sinbeta);
        y1 = y0 + (major * cosalpha * sinbeta + minor * sinalpha * cosbeta);
        plot(x1,y1,co,'LineWidth',lw);hold on; 
        plot(hex_peaks(:,1),hex_peaks(:,2),'ok','LineWidth',lw);hold on; %,'MarkerSize', 20
        xMajor1 = x0 + major * cos(phi); xMajor2 = x0 - major * cos(phi);
        yMajor1 = y0 + major * sin(phi); yMajor2 = y0 - major * sin(phi);
        p1=xMajor1:(xMajor2-xMajor1)/10:xMajor2; p2=yMajor1:(yMajor2-yMajor1)/10:yMajor2;
        if ~isempty(p1) & ~isempty(p2)
            plot(p1,p2, co,'LineWidth',lw); hold on;
        end
        xMinor1 = x0 + minor * cos(phi+pi/2); xMinor2 = x0 - minor * cos(phi+pi/2);
        yMinor1 = y0 + minor * sin(phi+pi/2); yMinor2 = y0 - minor * sin(phi+pi/2);
        p11=xMinor1:(xMinor2-xMinor1)/10:xMinor2; p21=yMinor1:(yMinor2-yMinor1)/10:yMinor2;
        if ~isempty(p11) & ~isempty(p21)
            plot(p11,p21, co,'LineWidth',lw); hold on;
            %REMOVE POINTS OUT OF RANGE<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        end
    end
    xlim(x); ylim(y);
end

function g = sortByMiddle(g)
   fs = fieldnames(g);
   a = struct2cell(g);
   s = size(a);
   a = reshape(a, s(1), []);
   a = a';
   for i = 1:length(a(:,1));
      a(i,11) = {length(a{i,9})}; %middle row
   end
   a = sortrows(a, -11);
   a = a(:,1:10);
   a = reshape(a', s);
   g = cell2struct(a, fs, 1);
end

function g = sortByEllipse(g)
   fs = fieldnames(g);
   a = struct2cell(g);
   s = size(a);
   a = reshape(a, s(1), []);
   a = a';
   for i = 1:length(a(:,1));
      %sort by ellipse size of *before*
      a(i,11) = {a{i,8}.('module').('major_ax') * a{i,8}.('module').('minor_ax')*pi};
      if isnan(a{i,8}.('module').('major_ax'))  || isnan(a{i,8}.('module').('minor_ax'))
          a(i,11) = {10000}; %make unknown max
      end
   end
   a = sortrows(a, 11); %-11 for backwards
   a = a(:,1:10);
   a = reshape(a', s);
   g = cell2struct(a, fs, 1);
end

function b=goodCell(c, rmt)
    b= c.exists && ...
    not(isnan(c.gridscore) || ...
         c.gridscore == -2 || ... 
         c.max_r < rmt     || ...
         c.max_r == 50        ...
     );
end

%{
remove corner cases
%}
function plot_group(g, key, iter, bs, m)
    prnt = true;
    rmt = 1.5; % Rate Map Threshold Hz
    %get dimensions of plot
    r = 6; %default plot length
    if (mod(length(g), r) ~= 0 & length(g) > 2*r) || length(g) == r+1
        r = r+1;
    end
    debug = '';
    %m = length(g(1).middle);
    r = 6; %m = 5; %m = # makes all plots same size 
    c = 2*(2+m) + 1; % width of figure: (before + max middle + after) * 2, for ratemap + autocorr + mcross
    %assums longest middle vector contains times of others
    %WHICH CELLS TO DISPLAY
    gridThresh =  0.7;
    good = [g(1)]; bad = [g(1)]; %find low gridscore cells
    for j = 1:length(g) %put low grid scores at end
        if g(j).before.gridscore > gridThresh
            good(end+1) = g(j);
        else
            bad(end+1) = g(j);
        end
    end
    g = good(2:end); %[good(2:end) bad(2:end)]; Only show good sells    
    if ~isempty(g) 
        g = sortByMiddle(g); %align times of cells
        %TODO: move to end if 0.3 gridness or less for before
        for j = 1:length(g(1).middle) %used to align middles
            t(j) = g(1).middle{j}.pt(1); %CHECK TIMES ALL ALIGN??
        end
        maxmid = min(length(t), m); %find index for after
        g = sortByEllipse(g); %display by ellipse size;
    end
    %PLOT
    tot = 0;
    page = 0; height = 0';
    while tot < length(g);
        page = page + 1;
        rr = r;                              %CHANGE r
        if r > (length(g) - tot)*2;  %*2
            rr = (length(g) - tot)*2; %*2
            height = round(300/r);
        end
        h = figure('Position', [0, 0, 2000, 1000]);%(2*m+5)*120 , 140*r + height]); %
        set(gca,'LooseInset', get(gca,'TightInset'));
        colormap jet;
        titl = sprintf('%dmin_i%d_%s_cells%d_pg%dof%d',bs,iter,key,length(g),page,ceil(2*length(g)/r));
        set(suptitle(titl),'Interpreter', 'none'); %PUT SUP TITLE AFTER ALL SUBPLOT COMMANDS << test
        %%% plot each cell %%%
        z = 0;
        while z < rr%rr %for
            z = z+1;
            tot = tot + 1;
            r1 = g(tot);
            mcross = {};
            %%%plot rm
            subplot(r, c, c*(z-1) + 1);
            gc = '';
            if r1.before.gridscore < gridThresh
                gc = '*';
            end
            imagesc(r1.before.rm), title(sprintf('%sc%d: t(-1) %.0fHz',gc, r1.ind, r1.before.max_r)); 
            mcross{1} = r1.before.rm; ind = 1 + m;
            %plot middle rms
            for i = 1: min(length(r1.middle), m);
                if  goodCell(r1.middle{i}, rmt)
                    %not(r1.middle{i}.max_r < rmt || r1.middle{i}.max_r == 50)
                    ind = find(t == r1.middle{i}.pt(1)); %align middles;
                    subplot(r, c, c*(z-1) + 1 + ind);
                    imagesc(r1.middle{i}.rm);
                    title(sprintf('%.0f-%.0f(%0.fHz)', r1.middle{i}.pt(1)/60,...
                        r1.middle{i}.pt(length(r1.middle{i}.pt))/60, r1.middle{i}.max_r));
                    %title(sprintf('%.0fHz', r1.middle{i}.max_r));
                    mcross{end + 1} = r1.middle{i}.rm;
                else
                 ind = find(t == r1.middle{i}.pt(1)); %align middles;
                 subplot(r, c, c*(z-1) + 1 + ind);
                 %imagesc(-1);
                 %title(sprintf('%.0f-%.0f(>%0.fHz)', r1.middle{i}.pt(1)/60,...
                        %r1.middle{i}.pt(length(r1.middle{i}.pt))/60, rmt));
                end
                
            end                           %after
            if r1.after.exists == 1 %will be fixed
                %subplot(r, c, c*(z-1) + 1 + m + 1); %always plot at fixed 
                subplot(r, c, c*(z-1) + 1 + maxmid + 1); %plot after last mid
                imagesc(r1.after.rm);title(sprintf('t(+1) %.0fHz', r1.after.max_r));
                mcross{end + 1} = r1.after.rm; 
            end
            %%%plot ac 
            plotAcorrModule(r1.before, h, [r, c, c*(z-1) + 3 + m], 't(-1):');
            %plot middle ac
            for i = 1: min(length(r1.middle), m);
                if goodCell(r1.middle{i}, rmt)
                    ind = find(t == r1.middle{i}.pt(1));
                    %subplot(r, c, c*(z-1) + 1 + ind);
                    plotAcorrModule(r1.middle{i}, h, [r, c, c*(z-1) + 3 + m + ind], '');
                else
                 ind = find(t == r1.middle{i}.pt(1)); %align middles;
                 subplot(r, c, c*(z-1) + 3 + m + ind); %(m,n,l);
                 %title(sprintf('%s%.2f(%d)','', r1.middle{i}.gridscore, length((r1.middle{i}.sx))));
                end
            end                           %after
            if r1.after.exists == 1  %will be fixed
                %subplot(r, c, c*(z-1) + 3 + 2*m + 1);
                %plotAcorrModule(r1.after, h, [r, c, c*(z-1) + 2 + 2*m + 2], 't(+1):'); %always plot last
                plotAcorrModule(r1.after, h, [r, c, c*(z-1) + 3 + m + maxmid + 1], 't(+1):'); %plot after mid 
                
            end
            %mcross
            subplot(r, c, c*(z-1) + 2 + 2*m + 3); mc = allCorr(mcross,mcross);
            mc =  padarray(mc,m+2-size(mc),-1,'post'); %HEREEHHE
            imagesc(mc, [0, 1]);  caxis([-1 1]);
            title(sprintf('[%.1f %.1f]', min(mc(:)), max(mc(:))));
            
           
            
            %%%%%%%%%%%%%%% TRAJECTORY
            z = z+1;
            subplot(r, c, c*(z-1) + 1);
            plot(r1.before.px, flip(r1.before.py));
            hold on
            scatter(r1.before.sx, flip(r1.before.sy), '.'),...
            xlim([0 100]), ylim([0 100]);
            title(sprintf('%.2f', r1.before.gridscore));
            mcrossT{1} = r1.before.rm;
            %plot middle rms
            for i = 1: min(length(r1.middle), m);
                if  goodCell(r1.middle{i}, rmt)
                    %not(r1.middle{i}.max_r < rmt || r1.middle{i}.max_r == 50)
                    ind = find(t == r1.middle{i}.pt(1)); %align middles;
                    subplot(r, c, c*(z-1) + 1 + ind);
                    plot(r1.middle{i}.px, flip(r1.middle{i}.py));
                    hold on
                    scatter(r1.middle{i}.sx, flip(r1.middle{i}.sy), '.'),...
                    xlim([0 100]), ylim([0 100]);
                    title(sprintf('%.2f', r1.middle{i}.gridscore));
                    %title(sprintf('%.0fHz', r1.middle{i}.max_r));
                    mcrossT{end + 1} = r1.middle{i}.rm;
                else
                     ind = find(t == r1.middle{i}.pt(1)); %align middles;
                     subplot(r, c, c*(z-1) + 1 + ind);
                end
                
            end                           %after
            if r1.after.exists == 1 %will be fixed
                %subplot(r, c, c*(z-1) + 1 + m + 1);
                subplot(r, c, c*(z-1) + 1 + maxmid + 1); %plot after last mid
                plot(r1.after.px, flip(r1.after.py));
                hold on
                scatter(r1.after.sx, flip(r1.after.sy), '.'),...
                xlim([0 100]), ylim([0 100]);
                title(sprintf('%.2f', r1.after.gridscore));
                mcrossT{end + 1} = r1.after.rm; 
            end

            %time correlation
 %           subplot(r, c, c*(z-1) + 3 + m);
 %           t2 = g(j).before.st;
 %           [p, c] = time_correlation(t1,t2);
            
            
            
            %MAKE AXIS NICE
             ax = findobj(gcf,'Type','Axes');
            for i=1:length(ax)
                set(ax(i),'FontSize',9);
                axis(ax(i),'equal')
                axis(ax(i),'off')
                %axis off; axis equal;
                %title(ax(i),{'Very Nice'})
            end
            
         
        end %% END ROW %%
        
        
        
        % END PLOTTING PAGE %%
        %[maxc maxl] = plotGroupTemporralCorrs(g, 500);
        
        %subplot(r,c,r*c);
        %plot(1,1);
        %axis off; axis equal;
        % ********* 
        sdir = 'C:\Noam\Output\muscimol\groups\';
        filename = sprintf('%s%s%s.png',sdir, debug, titl);
        %titl = sprintf('cells%d_pg%dof%d\n', length(g), page, ceil(length(g)/r));
        disp(filename);
        h.PaperPositionMode = 'auto'; %fig %fig = gcf;
        if prnt 
            print(filename, '-dpng','-r0');
            close(h);
        end
        %}
    end %END PLOTTING GROUP
    %}
end


function groups = find_simultaneously_recorded_cells(cells)
    groups = [];
    for i=1:length(cells)
        a = cells{i};
        if ~isempty(a.middle) %skip cells with non middle
            key = sprintf('g_%s_%s', a.id, a.date);
            if isfield(groups, key)
                t = groups.(key);
                t = [t a];
                groups.(key) = t;
            else
                groups = setfield(groups, key, a);
            end
        end
        %    fprintf('%d: bad cell %f rmthr %f', i, a.max_r, rate_map_thresh);
    end

    %validate
    keys = fieldnames(groups);
    for i = 1:length(keys)
        k = keys{i};
        group = groups.(k);
        t = group(1);
        t = length(t.before.px);
        for j = 1:length(group) %CHECK TIMES !!
            tt = group(j); %cells{group(j)};
            tt = length(tt.before.px);
            if t ~= tt
                fprintf('%f %f whats up with this group b.length %s?\n',t, tt, k);
            end
        end
    end

    fprintf('grouped\n');
end

function plot_pair(r1, r2, id)
    figure('Position', [100, 100, 1200, 900]);
    set(gca,'LooseInset',get(gca,'TightInset'));
    colormap jet; bad = '';
    rmt = 0.1; % Rate Map Threshold
    %%% set(gca,'YDir','normal') instead of flip???
    
    %n = min(length(r2.middle), length(r1.middle));
    if length(r2.middle) ~= length(r1.middle)
        %fprintf('%d: error plotting, middle sessions different lengths', id);%bad = '*';
    end

    set(suptitle(sprintf('%sP%d(i%d,i%d): rat %s %s C1(t%dc%d %s %s) C2(t%dc%d %s %s)',...
        bad,id,r1.ind,r2.ind,r1.id,r1.date,r1.tet,r1.cel,r1.a,r1.type,r2.tet,r2.cel,r2.a,r2.type)),...
        'Interpreter', 'none');

    %line up bin times
    i = 1; i1 = 1; i2 = 1; mid1 = {}; mid2 = {};
    while i1 <= length(r1.middle) && i2 <= length(r2.middle)
        if r1.middle{i1}.pt(1) == r2.middle{i2}.pt(1) %Will only show traj,ratemap/ac for good ratemaps
            if not(r1.middle{i1}.max_r < rmt || r1.middle{i1}.max_r == 50 ...
                    || r2.middle{i2}.max_r < rmt || r2.middle{i2}.max_r == 50)
                mid1{i} = r1.middle{i1};
                mid2{i} = r2.middle{i2};
                i = i + 1;
            end
            i1 = i1 + 1;
            i2 = i2 + 1;
        elseif r1.middle{i1}.pt(1) < r2.middle{i2}.pt(1)
            i1 = i1 + 1;
        else
            i2 = i2 + 1;
        end
    end
    r1.middle = mid1;
    r2.middle = mid2;
    m = length(mid1);
    n = m +1+1; %how many cols in figure, +1 before +1 mcross
    after = 0; %display after;
    if not(r1.after.max_r < rmt || r1.after.max_r == 50 ...
            ||r2.after.max_r < rmt || r2.after.max_r == 50)
        n = n + 1;
        after = 1;
    end
    r = 8; %Rows for figure (trajectory ratemap autoc xcross axcross)


    %%%plot trajectory
    subplot(r, n, 1 + n*0);
    plot(r1.before.px, flip(r1.before.py));
    hold on
    scatter(r1.before.sx, flip(r1.before.sy), '.'),...
        xlim([0 100]), ylim([0 100]), ...
        title(sprintf('c1: before %d',length(r1.before.st)));
    axis off; axis equal;
    subplot(r, n, 1 + n*1);
    plot(r2.before.px, flip(r2.before.py));
    hold on
    scatter(r2.before.sx, flip(r2.before.sy), '.'),...
        xlim([0 100]), ylim([0 100]), axis off; axis equal;
    title(sprintf('c2: %d',length(r2.before.st))); %TITLE
    for i = 1:m;
        %if length(r1.middle{i}.px) > 1 && length(r2.middle{i}.px) > 1 %Display if more than one spike
        subplot(r, n, i+1 + n*0);
        plot(r1.middle{i}.px, flip(r1.middle{i}.py));
        hold on
        scatter(r1.middle{i}.sx, flip(r1.middle{i}.sy), '.'); xlim([0 100]), ylim([0 100]);
        %title(sprintf('%.0f-%.0f(%d)', min(r1.middle{i}.pt)/60,...
        %max(r1.middle{i}.pt)/60, length(r1.middle{i}.st)));
        title(sprintf('%.0f-%.0f(%d)', r1.middle{i}.pt(1)/60,...
            r1.middle{i}.pt(length(r1.middle{i}.pt))/60, length(r1.middle{i}.st)));
        axis off; axis equal;
        subplot(r, n, (i+1) + n*1);
        plot(r2.middle{i}.px, flip(r2.middle{i}.py));
        hold on
        scatter(r2.middle{i}.sx, flip(r2.middle{i}.sy), '.'), xlim([0 100]), ylim([0 100]);
        %title(sprintf('%d',length(r2.middle{i}.st)));
        title(sprintf('%.0f-%.0f(%d)', r2.middle{i}.pt(1)/60,...
            r2.middle{i}.pt(length(r2.middle{i}.pt))/60, length(r2.middle{i}.st)));
        axis off; axis equal;
    end
    %if length(r1.after.px) > 1 && length(r2.after.px) > 1
    if after
        subplot(r, n, n*1-1);
        plot(r1.after.px, flip(r1.after.py));
        hold on
        scatter(r1.after.sx, flip(r1.after.sy), '.'),...
            xlim([0 100]), ylim([0 100]), axis off; axis equal;
        title(sprintf('after %d',length(r1.after.st)));
        subplot(r, n, n*2-1);
        plot(r2.after.px, flip(r2.after.py));
        hold on
        scatter(r2.after.sx, flip(r2.after.sy), '.'),...
            xlim([0 100]), ylim([0 100]); axis off; axis equal;
        title(sprintf('%d',length(r2.after.st)));
    end

    %%%plot rm
    subplot(r, n, 1 + n*2);
    imagesc(r1.before.rm), title(sprintf('c1: %.0fHz', r1.before.max_r)); axis off; axis equal;
    subplot(r, n, 1 + n*3);
    imagesc(r2.before.rm), title(sprintf('c2: %.0fHz', r2.before.max_r)); axis off; axis equal;
    for i = 1: m;
        subplot(r, n, i+1 + n*2);
        imagesc(r1.middle{i}.rm);title(sprintf('%.0fHz', r1.middle{i}.max_r));axis off; axis equal;
        subplot(r, n, i+1 + n*3);
        imagesc(r2.middle{i}.rm);title(sprintf('%.0fHz', r2.middle{i}.max_r));axis off; axis equal;
    end                           %after
    if after
        subplot(r, n, n*3-1);
        imagesc(r1.after.rm);title(sprintf('%.0fHz', r1.after.max_r));axis off; axis equal;
        subplot(r, n, n*4-1);
        imagesc(r2.after.rm);title(sprintf('%.0fHz', r2.after.max_r));axis off; axis equal;
    end
    % mcross
    t1{1} = r1.before.rm; t2{1} = r2.before.rm;
    for i = 1:length(mid1)
        t1{i+1} = mid1{i}.rm; t2{i+1} = mid2{i}.rm;
    end
    if after
        t1{m+2} = r1.after.rm; t2{m+2} = r2.after.rm;
    end
    subplot(r, n, n*3);
    imagesc(allCorr(t1,t1));title('cross per');axis off; axis equal;
    subplot(r, n, n*4);
    imagesc(allCorr(t2,t2));title('cross per');axis off; axis equal;

    %%%plot ac
    subplot(r, n, 1 + n*4);
    imagesc(r1.before.ac), title(sprintf('c1: %.2f', r1.before.gridscore));
    axis off; axis equal;
    subplot(r, n, 1 + n*5);
    imagesc(r2.before.ac), title(sprintf('c2: %.2f', r2.before.gridscore));
    axis off; axis equal;
    for i = 1: m;
        subplot(r, n, i+1 + n*4);
        imagesc(r1.middle{i}.ac); title(sprintf('%.2f', r1.middle{i}.gridscore)); axis off; axis equal;
        subplot(r, n, i+1 + n*5);
        imagesc(r2.middle{i}.ac); title(sprintf('%.2f', r2.middle{i}.gridscore)); axis off; axis equal;
    end
    if after
        subplot(r, n, n*5-1);
        imagesc(r1.after.ac); axis off; axis equal;
        title(sprintf('%.2f', r1.after.gridscore));
        subplot(r, n, n*6-1);
        imagesc(r2.after.ac); axis off; axis equal;
        title(sprintf('%.2f', r2.after.gridscore));
    end
    % mcross
    t1 = {}; t2 = {}; t1{1} = r1.before.ac; t2{1} = r2.before.ac;
    for i = 1:length(mid1)
        t1{i+1} = mid1{i}.ac; t2{i+1} = mid2{i}.ac;
    end
    if after
        t1{m+2} = r1.after.ac; t2{m+2} = r2.after.ac;
    end
    subplot(r, n, n*5);
    imagesc(allCorr(t1,t1)); title('cross per');axis off; axis equal;
    subplot(r, n, n*6);
    imagesc(allCorr(t2,t2));title('cross per');axis off; axis equal;

    %%% cross corr
    subplot(r, n, 1 + n*6);
    cc = Cross_Correlation(r1.before.rm, r2.before.rm);
    imagesc(cc); title(sprintf('Xcorr: %s',''));axis off; axis equal;
    %acc
    subplot(r, n, 1 + n*7);
    acc = Cross_Correlation(cc, cc); acc = getSubmatrixFromCenter(acc, cc);
    imagesc(acc);
    title(sprintf('AcXcorr: %.2f', gridscore(acc, -1))); axis off; axis equal;
    for i = 1: m;
        subplot(r, n, i+1 + n*6);
        cc = Cross_Correlation(r1.middle{i}.rm, r2.middle{i}.rm);
        imagesc(cc);axis off; axis equal;
        %acc
        subplot(r, n, i+1 + n*7);
        acc = Cross_Correlation(cc, cc); acc = getSubmatrixFromCenter(acc, cc);
        imagesc(acc);axis off; axis equal;title(sprintf('%.2f', gridscore(acc, -1)));
    end

    if after
        subplot(r, n, n*7-1);
        cc = Cross_Correlation(r1.after.rm, r2.after.rm); imagesc(cc);axis off; axis equal;
        %acc
        subplot(r, n, n*8-1); acc = Cross_Correlation(cc, cc);
        getSubmatrixFromCenter(acc, cc); imagesc(acc);
        title(sprintf('%.2f', gridscore(acc, -1)));axis off; axis equal;
    end




    sdir = 'C:\Noam\Output\muscimol\pairs\';
    filename = sprintf('%s%d_Rat_%s_date_%s_C%d_t%d_c%d_C%d_t%d_c%d.png',sdir, id, r1.id, r1.date, r1.ind, r1.tet, r1.cel, r2.ind, r2.tet, r2.cel);
    disp(filename);
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(filename, '-dpng','-r0');
    close;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end plot
function f = flip(x)
    m = max(x);
    f = m-x;
end

%not perfect, submatrix will be a dimension higher than B, if odd
%(could fix using length instead of round(B)
function S = getSubmatrixFromCenter(A, B) %A big matrix, B little matrix
raca = round(size(A)/2);
rbcb = round(size(B)/2);
S = A((raca(1)-rbcb(1)+1:raca(1)+rbcb(1)),((raca(2)-rbcb(2)+1:raca(2)+rbcb(2))));
end

function A = allCorr(a, b) %
    A = zeros(length(a));
    for i = 1:length(a);
        for j = 1:length(b)
            t1 = a{i}(:); t2 = b{j}(:);
            %all this work to make two non-nan same length vectors!
            t1(isnan(t1))=0; t2(isnan(t2))=0;
            if length(t1) ~= length(t2)
                if mod(length(t1) + length(t2), 2) ~= 0
                    t2 = [t2 ; 0];
                end
                t1 = padarray(t1, [ceil((max([length(t1) length(t2)]) - (length(t1)))/2) 0]);
                t2 = padarray(t2, [ceil((max([length(t1) length(t2)]) - (length(t2)))/2) 0]);
            end
            %{
            if length(t1) > length(t2)
                %could use image stretch <<
                t2 = [t2; zeros(length(t1)-length(t2),1)];
            else
                t1 = [t1; zeros(length(t2)-length(t1),1)];
            end
            %}
            A(i,j) = corr(t1, t2); 
            if length(a{i}(:)) == length(b{j}(:))
                %A(i,j) = corr(a{i}(:), b{j}(:));
            else
                %A(i,j) = 1;
                %fprintf('wrong lengths %d %d\n',length(a{i}(:)),length(b{j}(:)));
            end
        end
    end
end



function r = process(data, binsize)
r.id   = data.rat;
r.date = data.date;
r.tet  = data.tetrode;
r.cell  = data.cell;
r.type = data.cell_type_obj_new;
r.a = data.area;
r.ind = data.ind;
p = data.B(1).pos_data;
s = data.B(1).spike_data;
%before
tr = rat_trajectory(p.x1, p.x2, p.y1, p.y2, p.t, s.x, s.x2, s.y, s.y2, s.ts);
r.before.px = tr.px; r.before.py = tr.py; r.before.pt = tr.pt;
r.before.sx = tr.sx; r.before.sy = tr.sy; r.before.st = tr.st;  
r.before.rm = Create_Rate_Map(tr.px, tr.py, tr.pt, tr.sx, tr.sy, tr.st);
r.before.ac = Cross_Correlation(r.before.rm, r.before.rm);
r.before.max_r = max(r.before.rm(:));
r.before.gridscore = gridscore(r.before.ac, r.ind);
r.before.module = Find_Module(r.before.ac);
r.before.exists = false;
if r.before.max_r ~= 0 && r.before.max_r ~= 50 && r.before.gridscore ~= -2
    r.before.exists = true;
end

%middle
p = data.B(2).pos_data;
s = data.B(2).spike_data;
tr = rat_trajectory(p.x1, p.x2, p.y1, p.y2, p.t, s.x, s.x2, s.y, s.y2, s.ts);
binLength = binsize; %minutes,
maxTime =  45;% minutes
bins = bin_trajactory(tr, binLength*60, maxTime*60); %BIN LENGTH
for i = 1:length(bins)
    t = bins{i};
    bins{i}.rm = Create_Rate_Map(t.px, t.py, t.pt, t.sx, t.sy, t.st);
    bins{i}.ac = Cross_Correlation(bins{i}.rm, bins{i}.rm);
    bins{i}.max_r = max(bins{i}.rm(:));
    %fprintf('%d %d\n', r.ind, i);
    bins{i}.gridscore = gridscore(bins{i}.ac, r.ind);
    bins{i}.module = Find_Module(bins{i}.ac);
    bins{i}.exists = false;
    if bins{i}.max_r ~= 0 && bins{i}.max_r ~= 50 && bins{i}.gridscore ~= -2
        bins{i}.exists = true;
    end
end
r.middle = bins;

%after
if length(data.B) == 3               %after exists
    p = data.B(3).pos_data;
    s = data.B(3).spike_data;
    tr = rat_trajectory(p.x1, p.x2, p.y1, p.y2, p.t, s.x, s.x2, s.y, s.y2, s.ts);
    r.after.px = tr.px; r.after.py = tr.py; r.after.pt = tr.pt;
    r.after.sx = tr.sx; r.after.sy = tr.sy; r.after.st = tr.st;  
    r.after.rm = Create_Rate_Map(tr.px, tr.py, tr.pt, tr.sx, tr.sy, tr.st);
    r.after.ac = Cross_Correlation(r.after.rm, r.after.rm);
    r.after.max_r = max(r.after.rm(:));
    r.after.gridscore = gridscore(r.after.ac, r.ind);
    r.after.exists = false;
    r.after.module = Find_Module(r.after.ac);
    if r.after.max_r ~= 0 && r.after.max_r ~= 50 && r.after.gridscore ~= -2
        r.after.exists = true;
    end
else                                %after doesn't exist
    tr = rat_trajectory([0],[0],[0],[0],[0],[0],[0],[0],[0],[0]);
    r.after.px = tr.px; r.after.py = tr.py; r.after.pt = tr.pt;
    r.after.sx = tr.sx; r.after.sy = tr.sy; r.after.st = tr.st;  
    r.after.rm = false; %Create_Rate_Map(tr.px, tr.py, tr.pt, tr.sx, tr.sy, tr.st);
    r.after.ac = false; %Cross_Correlation(r.after.rm, r.after.rm);
    r.after.max_r = -1; %max(r.after.rm(:));
    r.after.gridscore = -2; %gridscore(r.after.ac, r.ind);
    r.after.exists = false;
    r.after.module = false; 
end

end

function score = gridscore(ac, ind)
    try
        score = GridnessRadius(ac, FindROuter(ac));
    catch ME
        fprintf('%d: error calculting grid score %s\n', ind, ME.message);
        score = -2;
    end
end

function test()
    sa=0, sb= 0;
    for i = 1:length(bins)
        fprintf('%d: min %.3f max %.3f \n',i, min(bins{i}.pt)/60, max(bins{i}.pt)/60);
        sa = sa + length(bins{i}.pt);
        sb = sb + length(bins{i}.st);
    end
    fprintf('%d  %d \n',sa, sb);
end

%bin size in seconds, max time
%tr = trajectory; bs = bin size in seconds; mt = max time
function bins = bin_trajactory(tr, bs, mt)%
minSpikes = 10;%min  spikes to be a bin
minBinLength = 1; %in minutes
bins{1,1} = []; b = 1;
ind = 0; %index in current bin
t = tr.pt; st = tr.st;
%t0 = ceil(tr.pt(1)/(5*60))*5*60 - bs; %statrs at 5min intervals %tr.pt(1);
t0 = tr.pt(1);
for i = 1: length(t) %index in entire series
    if t(i) <= mt
        %breaks up bins based on bin length or 15min interval cutoffs
        if t(i) < t0 + bs %&& mod(t(i), bs) %mod breaks up bins on even intervals, non 0 means continue
            ind = ind +1;
            bins{b}.px(ind) = tr.px(i);
            bins{b}.py(ind) = tr.py(i);
            bins{b}.pt(ind) = tr.pt(i);
        else
            ind = 1;
            t0 = t(i);
            b = b + 1;
            bins{b}.px(ind) = tr.px(i);
            bins{b}.py(ind) = tr.py(i);
            bins{b}.pt(ind) = tr.pt(i);
        end
    end
end
%spikes
b = 1; ind = 0;
for i = 1:length(st)
    if st(i) <= mt && st(i) <= max(bins{length(bins)}.pt) %cuts off spikes past max trajectory time
        while st(i) > max(bins{b}.pt) %skips past bins with no spikes
            b = b + 1;
            ind = 0;
        end
        ind = ind +1;
        bins{b}.sx(ind) = tr.sx(i);
        bins{b}.sy(ind) = tr.sy(i);
        bins{b}.st(ind) = tr.st(i);
    end
end
%remove weak bins less than minspikes, min lengths 
b = 1;
t = {};
for i = 1:length(bins)
    %conditions
    if isfield(bins{i},'sx') == 1 && length(bins{i}.st) > minSpikes &&... %if isfield(bins{i},'sx') == 0
       floor((bins{i}.pt(end) - bins{i}.pt(1))/60) >= minBinLength %make bins at least this long
        %bins{i}.sx(1) = 0;
        %bins{i}.sy(1) = 0;
        %bins{i}.st(1) = min(bins{b}.pt);
        t{b} = bins{i};
        b = b + 1;
    end
end
bins = t;
end

function c = rat_trajectory(px1, px2, py1, py2, pt, sx1, sx2, sy1, sy2, st)
    c.px =  mean([px1, px2], 2)';
    c.px = c.px - min(c.px) + 0.00001; %no zeros
    c.py =  mean([py1, py2], 2)';
    c.py = c.py - min(c.py) + 0.00001;
    c.pt = pt;
    c.sx =  mean([sx1, sx2], 2)';
    c.sx = c.sx - min(c.sx) + 0.00001;
    c.sy =  mean([sy2, sy2], 2)';
    c.sy =c.sy - min(c.sy) + 0.00001;
    c.st = st;
    %{figure();plot(c.px, c.py);hold on; plot(c.sx, c.sy,'.');%}
end

function rate_mat = Create_Rate_Map(px, py, pt, sx, sy, st)
    posx = px; posy = py; post = pt;
    spkx = sx; spky = sy; spkt = st;
    parms.sigma = 3; % gaussiam smoothing factor
    parms.time_per_bin=0.02; %defult size of bin time sampeling (1/sampeling rate)(updated in Read_Examples_2 function)
    % Minimum radius used in the auto-correlogram when finding the best
    parms.bin_size = 3; % size of spacial bin (for create the rate map)
    parms.num_of_direction_bins=60; % for head-direction calculatino
    parms.max_lag=500; % max lag (in msec) for temporal autocorrelation
    parms.bin_num=50;

    min_x = min(floor(posx)); max_x = ceil(max(posx));  
    min_y = min(floor(posy)); max_y = ceil(max(posy));

    % divide the environment into spatial bins
    axis_x = min_x : parms.bin_size : max_x;
    axis_y = min_y : parms.bin_size : max_y;

    time_mat =  zeros(length(axis_y),length(axis_x));
    spike_mat = zeros(length(axis_y),length(axis_x));
    rate_mat =  zeros(length(axis_y),length(axis_x));

    %create time mat (compute how much time the rat spends in each bin)
    % find in each moment(time_per_bin) what spatial bin the rat is at and add the time_per_bin to
    for i = 1:length(post)
        if ~isnan(posx(i)) && ~isnan(posy(i))
            [min_val,x_ind] =  min(abs(posx(i) - axis_x));
            [min_val,y_ind] =  min(abs(posy(i) - axis_y));
            time_mat(y_ind,x_ind) = time_mat(y_ind,x_ind) + parms.time_per_bin;
            %if ~isnan(min_val)
        end
    end
    %create conut mat( count the num of spikes in each bin)
    for i = 1:length(spkt)
        if ~isnan(spkx(i)) && ~isnan(spky(i))
            [min_val,x_ind] =  min(abs(spkx(i) - axis_x));
            [min_val,y_ind] =  min(abs(spky(i) - axis_y));
            spike_mat(y_ind,x_ind)= spike_mat(y_ind,x_ind)+1;
            %if ~isnan(min_val)
        end
    end
    % create rate mat
    rate_mat=spike_mat./time_mat;
    rate_mat(rate_mat==inf)=NaN;
    rate_mat_unsmoothed = rate_mat;
    %create window
    h=fspecial('gaussian',3*[parms.sigma,parms.sigma],parms.sigma);
    rate_mat = nanconv2(rate_mat,h);
    rate_mat(isnan(rate_mat))=0; %bad??????nza
    disp('');
end

function pairs_by_file_index = find_pairs(cells)
pairs_by_file_index = [];
grid_cells = 0;
for i=1:length(cells) %a = load(strcat(params.dir_load,files(i).name));
    a = cells{i};
    %if (strcmp(a.cell_type_obj_new, 'GRID') || strcmp(a.cell_type_obj_new, 'CONJ')) % '_GRID' ?????
    if a.before.gridscore > 0.3
        for j=i+1:length(cells)
            b = cells{j};
            %adjust for missing third <<<<<<
            if (       b.before.gridscore > 0.3 ...
                    && strcmp(a.id, b.id) && strcmp(a.date, b.date)...
                    ... && (strcmp(a.cell_type_subj, 'GRID') || strcmp(a.cell_type_subj, 'CONJ'))  ...
                    && strcmp(a.a, b.a)  ...  %area????
                    && ( a.tet ~= b.tet || a.cell ~= b.cell) ...
                    ... && length(a.B) == 3 &&  length(b.B) == 3 ... %ignore cells without before/after
                    )
                %fprintf('same rat:%s date:%s a:%d,%d b:%d,%d \n', a.rat, ...
                %a.date, a.tetrode, a.cell, b.tetrode, b.cell);
                
                grid_cells = grid_cells + 1;
                pairs_by_file_index(grid_cells,1) = i;
                pairs_by_file_index(grid_cells,2) = j;
            end
        end
    else
        %disp (a.cell_type_subj);
    end
end
fprintf('%d pairs found\n',grid_cells);
end

function pairs_by_file_index = find_pairs_of_cells_mosimol(data)
pairs_by_file_index = [];
grid_cells = 0;
for i=1:length(data) %a = load(strcat(params.dir_load,files(i).name));
    a = data(i); a = a{1}.db;
    if (strcmp(a.cell_type_obj_new, 'GRID') || strcmp(a.cell_type_obj_new, 'CONJ')) % '_GRID' ?????
        for j=i+1:length(data)
            b = data(j); b = b{1}.db;
            if length(a.B) ~= 3
                disp('missing 3rd');
            end
            if ( strcmp(a.rat, b.rat) && strcmp(a.date, b.date)...
                    && (strcmp(a.cell_type_subj, 'GRID') || strcmp(a.cell_type_subj, 'CONJ'))  ...
                    && strcmp(a.area, b.area)  ...  %area????
                    && ( a.tetrode ~= b.tetrode || a.cell ~= b.cell) ...
                    && length(a.B) == 3 &&  length(b.B) == 3 ... %ignore cells without before/after
                    )
                %fprintf('same rat:%s date:%s a:%d,%d b:%d,%d \n', a.rat, ...
                %a.date, a.tetrode, a.cell, b.tetrode, b.cell);
                
                grid_cells = grid_cells + 1;
                pairs_by_file_index(grid_cells,1) = i;
                pairs_by_file_index(grid_cells,2) = j;
            end
        end
    else
        %disp (a.cell_type_subj);
    end
end
fprintf('%d paris found\n',grid_cells);
end

function R_outer = FindROuter(acorr)
if length(acorr) == 1 || isnan(max(acorr(:)));
    R_outer = -1;
    %return;   %EXIT
end
% calculate all the extrema points in the spatial autocorrelation
[zmax,imax,zmin,imin]= Extrema2(acorr);
[i,j]=ind2sub(size(acorr),imax);
%put all extrema points in dist
dist(:,1)=j;
dist(:,2)=i;
n=length(i);
%calculate the distance of all extrema to the central peak and put them in
%column 3
dist(1:n,3)=sqrt(  (i(1:n)-i(1)).^2 + (j(1:n)-j(1)).^2);
% sort the hexonal peaks by distance to the centeral peak
[score,ind]=sort(dist(:,3));
dist=dist(ind,:);
%zmax=zmax(ind);
R=dist(2,3);
count=1;
i=2;
hex_peaks(1,:,:)=dist(1,:,:);
% finds the first 6 closest peaks to the central peak
while count<7 && i<=length(dist)
    % calculate the min distance of point i from all other already chosen
    % points
    min_dist_peaks=min(sqrt(  (hex_peaks(1:size(hex_peaks,1),1)-dist(i,1)).^2 +...
        (hex_peaks(1:size(hex_peaks,1),2)-dist(i,2)) .^2));
    % point i needs to be on the right side (we choos only half the point cause its semetrical)
    % and the distance of point i from all other already chosen points
    % needs to be higher than R/2
    if dist(i,1)>=dist(1,1) && min_dist_peaks>(R/1.5)
        hex_peaks(count,:,:)=dist(i,:,:);
        count=count+1;
        hex_peaks(count,1)=dist(1,1)-dist(i,1)+dist(1,1);
        hex_peaks(count,2)=dist(1,2)-dist(i,2)+dist(1,2);
        
        count=count+1;
    end
    i=i+1;
end
R_outer=max(hex_peaks(:,3))*1.15;
end

function gridness2 = GridnessRadius(org_Acor ,R_outer)
if length(org_Acor) == 1 || isnan(max(org_Acor(:)))
    gridness2 = -1;
    %return; EXIT
end
% compute distance from center
[val,center_y]=max(max(org_Acor));
[val,center_x]=max(max(org_Acor'));
[Y,X] = ndgrid(1:1:size(org_Acor,1), 1:1:size(org_Acor,2));
dist_from_center=sqrt((Y-center_y).^2+(X-center_x).^2);
% making sure that outer radius of the ring (R_outer) is not bigger than the distance matrix
R_outer=min([min(dist_from_center(1,:)),min(dist_from_center(:,1)),...
    min(dist_from_center(size(dist_from_center,1),:))...
    min(dist_from_center(:,size(dist_from_center,2))),R_outer]);
% compute inner radius of the anulus (ring)
R_inner=ceil(min(dist_from_center(org_Acor<0.1)));
%extract the original anulus (ring) from Acor
org_Ring=org_Acor(dist_from_center<=R_outer & dist_from_center>=R_inner);
% make sure that after rotation and interpulation the center will remain the maximum point.
org_Acor(center_x,center_y)=10;
for jj = 2:6
    % rotate the auto-correlation
    rot_Acor=imrotate(org_Acor,(jj-1)*30,'bicubic');
    % compute distance from new center
    [val,tmp_center_x]=max(max(rot_Acor));
    [val,tmp_center_y]=max(max(rot_Acor'));
    [Y,X] = ndgrid(1:1:size(rot_Acor,1), 1:1:size(rot_Acor,2));
    tmp_dist_from_center=sqrt((Y-tmp_center_y).^2+(X-tmp_center_x).^2);
    % extract the anulus(ring)
    rot_Ring=rot_Acor(tmp_dist_from_center<=R_outer & tmp_dist_from_center>=R_inner);
    if length(rot_Ring)~=length(org_Ring)
        gridness2=nan;
        return
    end
    % compute pearson correlation between rotate Acor and original Acor
    corrValues(jj) = PointCorr(org_Ring,rot_Ring);
    clear rot_Ring tmp_center_x tmp_center_y tmp_dist_from_center Y X
end
% min of higher correlation at 60 and 120 degree rotation
min_rot_60_120 = min([corrValues(3),corrValues(5)]);
% max of lower correlations 30,90,150 rotation
max_rot_30_90_150 = max([corrValues(2),corrValues(4),corrValues(6)]);
% calculate gridness min(60,120)-max(30,90,150)
gridness2 = min_rot_60_120 - max_rot_30_90_150;
% different way to calculate gridness
gridness1=mean(([corrValues(3),corrValues(5)]))-...
    mean([corrValues(2),corrValues(4),corrValues(6)]);
end

function out_mat=PointCorr(Rxx,RxxR)
nan_mat=Rxx .* RxxR;
notnan_inds = find(~isnan(nan_mat));  %normalized to the number of nontnan components (average)
n=length(notnan_inds);
if n < 2
    out_mat = NaN;
end
sigma_x_y =sum(nan_mat(notnan_inds));
sigma_x =      sum(Rxx(notnan_inds));
sigma_y =      sum(RxxR(notnan_inds));
sigma_x2 = sum(Rxx(notnan_inds).^2);
sigma_y2 = sum(RxxR(notnan_inds).^2);

out_mat = (n*sigma_x_y - sigma_x.*sigma_y) ./ ...
    sqrt(n*sigma_x2-sigma_x.^2) ./ ...
    sqrt(n*sigma_y2-sigma_y.^2);
end




function out_mat=Cross_Correlation(mat1,mat2)

[ma,na] = size(mat1);
[mb,nb] = size(mat2);
%mc = max([ma+mb-1,ma,mb]);
%nc = max([na+nb-1,na,nb]);
%m=min(mc,nc);
%out_mat = nan(m,m);
out_mat = nan(ma+mb-1,na+nb-1);

i_size = size(mat2,1); j_size = size(mat2,2);
[work_mat,npad_i,npad_j] = pad_edges(mat1,mat2,1);

for i = 1:size(out_mat,1)
    for j = 1:size(out_mat,2)
        
        % for each i and j, choose the correct sub-mat (size of mat 2) to
        % multiply with mat2
        
        sub_mat = work_mat(npad_i+i-floor(i_size):npad_i+i-1, ...
            npad_j+j-floor(j_size):npad_j+j-1  );
        nan_sub_mat=sub_mat .* mat2;
        notnan_inds = find(~isnan(nan_sub_mat));  %normalized to the number of nontnan components (average)
        
        n=length(notnan_inds);
        
        if n < 20
            out_mat(i,j) = NaN;
            continue;
        end
        
        sigma_x_y =sum(nan_sub_mat(notnan_inds));
        sigma_x =      sum(sub_mat(notnan_inds));
        sigma_y =      sum(mat2(notnan_inds));
        sigma_x2 = sum(sub_mat(notnan_inds).^2);
        sigma_y2 = sum(mat2(notnan_inds).^2);
        
        out_mat(i,j) = (n*sigma_x_y - sigma_x.*sigma_y) ./ ...
            sqrt(n*sigma_x2-sigma_x.^2) ./ ...
            sqrt(n*sigma_y2-sigma_y.^2);
    end % for j
end % for i
disp('')
end




function out_mat = nanconv2(mat,h)
out_mat = mat;
nan_mat = isnan(mat);
% dilate nan_mat
SE = strel('disk', 2);
nan_mat =  ~imdilate(~nan_mat,SE);
i_size = size(h,1); j_size = size(h,2);
[work_mat,npad_i,npad_j] = pad_edges(mat,h,2);
for i = 1:size(mat,1)
    for j = 1:size(mat,2)
        % for each i and j, choose the correct sub-mat (size of h) to multiply with h
        sub_mat = work_mat(npad_i+i-floor(i_size/2):npad_i+i+floor(i_size/2), ...
            npad_j+j-floor(j_size/2):npad_j+j+floor(j_size/2)  ); % assumes h is odd in number
        
        notnan_inds = find(~isnan(sub_mat));
        if ~isempty(notnan_inds)
            sum_h = sum(h(notnan_inds));   % normalize to the places without a NaN
            out_mat(i,j) = nansum(nansum(sub_mat .* h));
            out_mat(i,j) = out_mat(i,j)/sum_h;
        end
    end % for j
end % for i
out_mat(nan_mat) = NaN;
disp('')
end


function [out_mat,npad_i,npad_j] = pad_edges(mat,h,l)
npad_ij = ceil(size(h)/l);
npad_i = npad_ij(1);
npad_j = npad_ij(2);
in_size = size(mat);
out_size = in_size + [2*npad_i 2*npad_j];
out_mat = nan(out_size);
out_mat(npad_i+1:npad_i+in_size(1),npad_j+1:npad_j+in_size(2)) = mat;
disp('')
end


function [semimajor_axis, semiminor_axis, x0, y0, phi] = ellipse_fit(x, y)
%
% Programmed by: Tal Hendel <thendel@tx.technion.ac.il>
% Faculty of Biomedical Engineering, Technion- Israel Institute of Technology
% 12-Dec-2008

x = x(:);
y = y(:);

%Construct M
M = [2*x.*y y.^2 2*x 2*y ones(size(x))];

% Multiply (-X.^2) by pseudoinverse(M)
e = M\(-x.^2);

%Extract parameters from vector e
a = 1;
b = e(1);
c = e(2);
d = e(3);
f = e(4);
g = e(5);

%Use Formulas from Mathworld to find semimajor_axis, semiminor_axis, x0, y0
%, and phi

delta = b^2-a*c;

x0 = (c*d - b*f)/delta;
y0 = (a*f - b*d)/delta;

phi = 0.5 * acot((c-a)/(2*b));

%phi = 0.5 * atan2((2*b),(c-a));

nom = 2 * (a*f^2 + c*d^2 + g*b^2 - 2*b*d*f - a*c*g);
s = sqrt(1 + (4*b^2)/(a-c)^2);

a_prime = sqrt(nom/(delta* ( (c-a)*s -(c+a))));

b_prime = sqrt(nom/(delta* ( (a-c)*s -(c+a))));

semimajor_axis = max(a_prime, b_prime);
semiminor_axis = min(a_prime, b_prime);

if (a_prime < b_prime)
    phi = pi/2 - phi;
end

disp('')
end


function module = Find_Module(acorr)

%REBEKKAHS
hex_peaks = find_six_points(acorr);

if hex_peaks ~= -1
    
    % fitting the elipse to the 6 peaks that we have found
    [major_ax, minor_ax, x0, y0, angle] = ellipse_fit(hex_peaks(1:6,1), hex_peaks(1:6,2));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% fixing a bug in ellipse fit %%%%%%%%%%%%%%%%
    [x1,y1,fig]=ellipse(major_ax,minor_ax,angle,x0,y0,'k',100);
    [x2,y2,fig]=ellipse(major_ax,minor_ax,-angle,x0,y0,'k',100);
    
    %     %%%%%%%%%%%%%%%%%%%% calculate which ellipse fits better (distance to points)
    i=1:length(x1);
    d1=0;d2=0;
    for j=1:6
        d1=d1+min((hex_peaks(j,1)-x1(i)).^2+(hex_peaks(j,2)-y1(i)).^2);
        d2=d2+min((hex_peaks(j,1)-x2(i)).^2+(hex_peaks(j,2)-y2(i)).^2);
    end
    %%%%%%%%%%%%%%%%% if d2 smaller put -angel and the bug is fixed
    if d2<d1
        angle=-angle;
    end
    
    if ~isreal(major_ax) ||  ~isreal(minor_ax)
        major_ax=nan;
        minor_ax=nan;
    end
    
    hex_peaks(7,1)=x0;
    hex_peaks(7,2)=y0;
    
    module.major_ax = major_ax;
    module.minor_ax = minor_ax;
    module.angle = angle;
    module.hex_peaks = hex_peaks;
    module.x0 = x0;
    module.y0 = y0;
    module.exists = 1;
    
else %ERROR CONDITION
    module.major_ax = -1;
    module.minor_ax = -1;
    module.angle = -1;
    module.hex_peaks = -1;
    module.x0 = -1;
    module.y0 = -1;
    module.exists = 0;
end
disp('')
end

function return_value=Distance(x,y,x2,y2)
return_value = sqrt((x2-x).^2+(y2-y).^2);
end

function [six_orientation_pts] = find_six_points(autocorr)

%find autocorrelation mat maximum points
[auto_max_inds] = FindMaxIndsRateMap(autocorr);

if length(auto_max_inds) >= 7
    
    [size_x, size_y]=size(autocorr);
    
    %finds distances from center of auto_corr_map to all max peaks
    cen = 1:length(auto_max_inds);
    auto_distances = Distance(auto_max_inds(cen, 1),auto_max_inds(cen, 2),(size_x/2)+0.5,(size_y/2)+0.5);
    
    
    
    %find point closest to center
    min_distance_index = find(auto_distances==min(auto_distances));
    middle_pt(1) = auto_max_inds(min_distance_index(1),1); %NZA Added(1) because rarely there are two vals
    middle_pt(2) = auto_max_inds(min_distance_index(1),2); %same
    
    %repeats find distances with more accurate center point
    cen = 1:length(auto_max_inds);
    
    auto_distances = Distance(auto_max_inds(cen, 1),auto_max_inds(cen, 2),middle_pt(1),middle_pt(2));
    
    [new_distances, six_inds] = sort(auto_distances (:));
    
    % finds the points for the module ellipse
    % finds the 6 closest pts to the center
    auto_dist_inds1 = find(auto_distances == new_distances(2));
    auto_dist_inds2 = find(auto_distances == new_distances(3));
    auto_dist_inds3 = find(auto_distances == new_distances(4));
    auto_dist_inds4 = find(auto_distances == new_distances(5));
    auto_dist_inds5 = find(auto_distances == new_distances(6));
    auto_dist_inds6 = find(auto_distances == new_distances(7));
    
    union1 = union(auto_dist_inds1, auto_dist_inds2);
    union2 = union(auto_dist_inds3, auto_dist_inds4);
    union3 = union(auto_dist_inds5, auto_dist_inds6);
    union1 = union(union1, union2);
    auto_dist_inds = union(union1, union3);
    
    %
    six_orientation_pts=nan(length(auto_dist_inds),2);
    for k= 1:length(auto_dist_inds);
        six_orientation_pts (k,1) = auto_max_inds(auto_dist_inds(k),1);
        six_orientation_pts (k,2) = auto_max_inds(auto_dist_inds(k),2);
    end
    %}
    
    %adds center point
    
    six_orientation_pts (7,1) = middle_pt(1);
    six_orientation_pts (7,2) = middle_pt(2);
    
    % to check for accuracy:
    % figure; imagesc(autocorr); hold on;
    % plot(six_orientation_pts(:,2), six_orientation_pts(:,1), 'x')
    
    disp('');
else
    six_orientation_pts = -1;
end

end

function [max_inds] = FindMaxIndsRateMap(rate_mat)

max_inds = 0;

% turn nans in rate_mat to zeros
rate_mat(isnan(rate_mat))=0;

% pad rate mat with zero edges to catch border maxs
rm_size= size(rate_mat);
rm_size= rm_size+2;
rate_mat_new= zeros(rm_size);
rate_mat_new(2:end-1, 2:end-1)= rate_mat(1:end, 1:end);
rate_mat=rate_mat_new;

max_inds_len = 0;

[size_x, size_y]= size(rate_mat);

for fig_i = 2:size_x-1
    for j = 2:size_y-1
        if rate_mat(fig_i,j) > rate_mat(fig_i+1,j) && ...
                rate_mat(fig_i,j) > rate_mat(fig_i-1,j) && ...
                rate_mat(fig_i,j) > rate_mat(fig_i,j+1) && ...
                rate_mat(fig_i,j) > rate_mat(fig_i,j-1)
            hold on
            max_inds_len = max_inds_len+1;
            max_inds(max_inds_len,1) = fig_i-1;     %list of maximum pt indices
            max_inds(max_inds_len,2) = j-1;
        end
    end
end
end

function [x,y,h]=ellipse(ra,rb,ang,x0,y0,C,Nb)
% Ellipse adds ellipses to the current plot
%
% ELLIPSE(ra,rb,ang,x0,y0) adds an ellipse with semimajor axis of ra,
% a semimajor axis of radius rb, a semimajor axis of ang, centered at
% the point x0,y0.
%
% The length of ra, rb, and ang should be the same.
% If ra is a vector of length L and x0,y0 scalars, L ellipses
% are added at point x0,y0.
% If ra is a scalar and x0,y0 vectors of length M, M ellipse are with the same
% radii are added at the points x0,y0.
% If ra, x0, y0 are vectors of the same length L=M, M ellipses are added.
% If ra is a vector of length L and x0, y0 are  vectors of length
% M~=L, L*M ellipses are added, at each point x0,y0, L ellipses of radius ra.
%
% ELLIPSE(ra,rb,ang,x0,y0,C)
% adds ellipses of color C. C may be a string ('r','b',...) or the RGB value.
% If no color is specified, it makes automatic use of the colors specified by
% the axes ColorOrder property. For several circles C may be a vector.
%
% ELLIPSE(ra,rb,ang,x0,y0,C,Nb), Nb specifies the number of points
% used to draw the ellipse. The default value is 300. Nb may be used
% for each ellipse individually.
%
% h=ELLIPSE(...) returns the handles to the ellipses.
%
% as a sample of how ellipse works, the following produces a red ellipse
% tipped up at a 45 deg axis from the x axis
% ellipse(1,2,pi/8,1,1,'r')
%
% note that if ra=rb, ELLIPSE plots a circle
%

% written by D.G. Long, Brigham Young University, based on the
% CIRCLES.m original
% written by Peter Blattner, Institute of Microtechnology, University of
% Neuchatel, Switzerland, blattner@imt.unine.ch


% Check the number of input arguments

if nargin<1,
    ra=[];
end;
if nargin<2,
    rb=[];
end;
if nargin<3,
    ang=[];
end;

%if nargin==1,
%  error('Not enough arguments');
%end;

if nargin<5,
    x0=[];
    y0=[];
end;

if nargin<6,
    C=[];
end

if nargin<7,
    Nb=[];
end

% set up the default values

if isempty(ra),ra=1;end;
if isempty(rb),rb=1;end;
if isempty(ang),ang=0;end;
if isempty(x0),x0=0;end;
if isempty(y0),y0=0;end;
if isempty(Nb),Nb=300;end;
if isempty(C),C=get(gca,'colororder');end;

% work on the variable sizes

x0=x0(:);
y0=y0(:);
ra=ra(:);
rb=rb(:);
ang=ang(:);
Nb=Nb(:);

if isstr(C),C=C(:);end;

if length(ra)~=length(rb),
    error('length(ra)~=length(rb)');
end;
if length(x0)~=length(y0),
    error('length(x0)~=length(y0)');
end;

% how many inscribed elllipses are plotted

if length(ra)~=length(x0)
    maxk=length(ra)*length(x0);
else
    maxk=length(ra);
end;

% drawing loop

for k=1:maxk
    
    if length(x0)==1
        xpos=x0;
        ypos=y0;
        radm=ra(k);
        radn=rb(k);
        if length(ang)==1
            an=ang;
        else
            an=ang(k);
        end;
    elseif length(ra)==1
        xpos=x0(k);
        ypos=y0(k);
        radm=ra;
        radn=rb;
        an=ang;
    elseif length(x0)==length(ra)
        xpos=x0(k);
        ypos=y0(k);
        radm=ra(k);
        radn=rb(k);
        an=ang(k)
    else
        rada=ra(fix((k-1)/size(x0,1))+1);
        radb=rb(fix((k-1)/size(x0,1))+1);
        an=ang(fix((k-1)/size(x0,1))+1);
        xpos=x0(rem(k-1,size(x0,1))+1);
        ypos=y0(rem(k-1,size(y0,1))+1);
    end;
    
    co=cos(an);
    si=sin(an);
    the=linspace(0,2*pi,Nb(rem(k-1,size(Nb,1))+1,:)+1);
    x=radm*cos(the)*co-si*radn*sin(the)+xpos;
    y=radm*cos(the)*si+co*radn*sin(the)+ypos;
    h(k)=line(radm*cos(the)*co-si*radn*sin(the)+xpos,radm*cos(the)*si+co*radn*sin(the)+ypos);
    set(h(k),'color',C(rem(k-1,size(C,1))+1,:));
    close; %noam
    
end
end

function misc()
    datastats(pos_by_bin);
    %datastats(spk_bins)
    % hist(a,9);figure;hist(b,9);
    dir_bins_pos = nan(nbins,length(pos_md));
%    dir_bins_spk = nan(nbins,length(spk_md));
    dir_bins_c1 = cell(nbins,1);
    dir_bins_c2 = cell(nbins,1);
    scale = [0:1e-3:max(max(c1.st),max(c2.st))];
    %one timestep contains multiple spikes
    tst1 = round(c1.st*1000)/1000;
    tst2 = round(c2.st*1000)/1000; %
    
    edges = [-Inf, mean([c1.pt(2:end) c1.pt(1:end-1)],2)', +Inf];
    I = discretize(c1.st, edges); %index of closest PT time for ST
    t = abs(c1.st - c1.pt(I)) > dt;
    
    for i = 1:length(c1.st)
        if min(abs(c1.pt-tst2(i)))<dt
        end
    end
    
    for i = 1:length(tst2)
        if min(abs(c1.pt-tst2(i)))<dt %find closest value to firing time within threshold
            dir_bins_c1{pos_by_bin(i)}(end+1) = tst1(i);
        end
    end
    
    for i = 1:length(pos_by_bin)
        dir_bins_pos(pos_by_bin(i),i) = c1.pt(i);
    end
    
    for i = 1:length(spk1_by_bin)
        dir_bins_spk(spk1_by_bin(i),i) = c1.st(i);
    end
    %    dir_bins_c1(pos_by_bin(i),end+1) = min(abs(c1.st-c1.pt(i)))<dt;
    %    dir_bins_c2(pos_by_bin(i),end+1) = min(abs(c2.st-c1.pt(i)))<dt;
    %
    %CHECK IF ALL SPIKE TIMES CONTAINED IN POS TIMES
    %map both sets of spikes to times it was in that direction
    
    
    tol=dt
    th = 1.6;
    [ii,jj]=find(abs(A-th)<tol)
    size(ii)
    
    %               generating trains
    %scales times x 1000, +1 so non-zero
    spkt1=floor(spkt1*1000)+1;
    spkt2=floor(spkt2*1000)+1;
    max_time=max(max(spkt1),max(spkt2));
    %generating the spike train
    train=zeros(2,max_time);
    train(1,spkt1)=1; %array at these indices(time*1000) will be 1, else 0
    train(2,spkt2)=1;
    
        %pos_md = atan2(vy,vx) + pi;%get rid of 0, 360 overlap %THINK ABOUT THIS
    %spk_md = interp1(c1.pt, unwrap(pos_md'), c1.st); 
    %spk_md = wrapToPi(spk_md) + pi; %get rid of 0, 360 overlap
    %rad_bins = 0:2*pi/(nbins):2*pi; %removed bins-1 because ends(pi -pi) are redundant
    %ang_bins = rad2deg(rad_bins);
    %count = hist(spk_MD,ang_ax); %spikes in each direction, bin are centers
    %count(end) = count(1) + count(end); count(1) = count(end);% noam, because ends are half bins, -180, 180
    %time = hist(pos_MD,ang_ax)*dt; %time in each direction
    %time(end) = time(1) + time(end); time(1) = time(end);% noam, because ends are half bins, -180, 180
    %basically I need to make spike trains of different directions
    %discretize(data, edges);
%             tst1 = round(dir_bins_spk1{i}*1000)+1;
%         tst2 = round(dir_bins_spk2{i}*1000)+1; %  
%         min_spk_time = min([tst1;tst2]);
%         max_spk_time = max(max([tst1;tst2]));
%         train1=zeros(min_spk_time,max_spk_time);
%         train2=zeros(min_spk_time,max_spk_time);
%         train1(tst1)=1;
%         train2(tst2)=1;
end


