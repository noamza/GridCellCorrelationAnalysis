function [pearson_xcov, count_xcor] = compareByMovingDirection(c1, c2, params)
    nbins = params.number_degree_bins;
    if nbins == 1
        pearson_xcov = cell(1,3); count_xcor = cell(1,3);
        spkt1=floor(c1.st*1000)+1; spkt2=floor(c2.st*1000)+1; min_time=min(min(spkt1),min(spkt2));
        spkt1 = spkt1-min_time+1; spkt2 = spkt2-min_time+1; max_time=max(max(spkt1),max(spkt2));
        train1=zeros(1,max_time);train2 =zeros(1,max_time); train1(spkt1)=1; train2(spkt2)=1;
        pearson_xcov{1,1}=xcov(train1,train2,round(1000*params.lag_max_secs),'coef')';
        [b,a] = butter(6,0.002*1.7^params.sigma); %LOW PASS 0.15%6th order, fc/fs/2 determined empiracally 
        pearson_xcov{1,1} = filtfilt(b,a,pearson_xcov{1,1});
        count_xcor{1,1}=xcorr(train1,train2, round(1000*params.lag_max_secs))'; %how many times spike at same point
        time_scale=( (1:length(pearson_xcov{1,1})) - ((length(pearson_xcov{1,1})-1)/2) - 1)/1000;
        pearson_xcov{1,2} = time_scale; count_xcor{1,3}= time_scale;
        pearson_xcov{1,3} = [0 360]; count_xcor{1,3}   = [0 360];
        return
    end
    px = double(c1.px); py = double(c1.py); pt = double(c1.pt);
    dt = median(diff(pt)); dt = round(dt*1e4)/1e4; % round to 10th of millisec
    vox = zeros(length(pt),1);
    voy = zeros(length(pt),1);
    %figure;subplot(2,2,1); imagesc(c1.rm); title('c1');subplot(2,2,2); imagesc(c2.rm);title('c2');
    %subplot(2,2,3); imagesc(xcorr2(c1.rm, c2.rm)); title('c1xc2'); colormap jet;
    for i = 2:length(vox);
        vox(i)=(px(i)-px(i-1))/(pt(i)-pt(i-1));
        voy(i)=(py(i)-py(i-1))/(pt(i)-pt(i-1));
    end %SMOOTH VELOCITY %test by reconstructing
    [b2,b1] = butter(6,0.06);%0.06);%0.1 
    vx = filtfilt(b2,b1,vox); 
    vy = filtfilt(b2,b1,voy);
    %REMOVE VELOCITY OVER 3cm/s
    for i = 1:length(vx)
        if sqrt((vx(i)).^2 + (vy(i)).^2) <= 3
            vx(i) = nan; vy(i) = nan;
        end
    end
    %figure();plot(diff(c1.px)/dt);hold on;plot(vx)
    pos_md = atan2(vy,vx);%wrapTo360(rad2deg(atan2d(vy,vx));
    %CHECK IF DELETING A LOT
    c1.st = c1.st(c1.st>=min(pt)&c1.st<=max(pt));
    c2.st = c2.st(c2.st>=min(pt)&c2.st<=max(pt));
    spk1_md = interp1(pt, unwrap(pos_md'),c1.st);
    spk2_md = interp1(pt, unwrap(pos_md'),c2.st);
    %TO 360 
    pos_md = wrapTo360(rad2deg(pos_md)); %check??
    spk1_md = wrapTo360(rad2deg(spk1_md)); spk2_md = wrapTo360(rad2deg(spk2_md));
    %BIN BY DIRECTIONS
    bin_edges = (0:360/nbins:360) - 180/nbins; bin_edges(1) = 0; bin_edges(end +1) = 360;% CHECK
    pos_by_bin = discretize(pos_md, bin_edges); pos_by_bin(pos_by_bin==nbins+1) = 1; %conbines two last bins into one, to wrap pi
    %t = abs(diff(pos_by_bin)); %check percent jump in opposite directions bin
    %[sum(t==2) sum(t==0) sum(t==1) sum(t==3) sum(isnan(t))]/length(t)*100;
    %Checked by comparing to pos bins, seems reasonable
    spk1_by_bin = discretize(spk1_md, bin_edges);spk1_by_bin(spk1_by_bin==nbins+1) = 1; %last bin added to first
    spk2_by_bin = discretize(spk2_md, bin_edges);spk2_by_bin(spk2_by_bin==nbins+1) = 1; 

    %TIME INTERVALS
    time_in_bin = cell(nbins,1);
    i = 1;
    while i <= length(pt)
        cb = pos_by_bin(i);
        increment = true;
        if ~isnan(cb)
            start = pt(i);
            while pos_by_bin(i) == cb && i < length(pt) %will enter first time without changing i
                ends = pt(i);
                i = i+1;
                increment = false;
            end
            start = round(start/dt)*dt; ends = round(ends/dt)*dt; %fancy rounding to dt!
            time_in_bin{cb} = [time_in_bin{cb}; [start:dt:ends]'];
        end
        if increment
            i = i+1;
        end
    end
    
    %BUILD TRAINS
    trains_dir_only=cell(nbins,1);
    for i=1:nbins
        trains_dir_only{i} = zeros(2,length(time_in_bin{i}));
    end
    %merge spikes into time array at closest time within dt
    for i=1:length(spk1_by_bin)
        if ~isnan(spk1_by_bin(i)) %which direction spike in    spike time
            [dif, idx] = min(abs(time_in_bin{spk1_by_bin(i)} - c1.st(i)));
            if dif <= dt
                trains_dir_only{spk1_by_bin(i)}(1, idx) = ... %more than one spike per time unit
                    trains_dir_only{spk1_by_bin(i)}(1, idx) + 1; 
            else
              %  fprintf('%d UNSORTED 1 SPIKE %f\n',i,dif);
            end
        end
    end
    for i=1:length(spk2_by_bin)
        if ~isnan(spk2_by_bin(i)) %which direction spike in    spike time
            [dif, idx] = min(abs(time_in_bin{spk2_by_bin(i)} - c2.st(i)));
            if dif <= dt
                trains_dir_only{spk2_by_bin(i)}(2, idx) = trains_dir_only{spk2_by_bin(i)}(2, idx) + 1;
            else
               % fprintf('%d UNSORTED 2 SPIKE %f\n',i,dif);
            end
        end
    end
    binned_trains = trains_dir_only;
    %BIN BY # OF SECS
    bin_size = round(params.time_bin_secs/dt);
    %ALWAYS ZEROS AT END??
    if bin_size > 1
        binned_trains = cell(4,1);
        for i=1:nbins
            z = trains_dir_only{i}';
            if mod(length(z),bin_size) ~= 0
                z = [z;zeros(bin_size-mod(length(z),bin_size),2)];
            end
            %reshapes to rectangle width of bin size, sums cols;
            binned_trains{i}(1,:) = sum(reshape(z(:,1),bin_size, ceil(length(z)/bin_size)));
            binned_trains{i}(2,:) = sum(reshape(z(:,2),bin_size, ceil(length(z)/bin_size))); 
        end
    end
    
    %CORRELATION
    pearson_xcov = cell(nbins,3); count_xcor = cell(nbins,3);
    lag_units = round(params.lag_max_secs/params.time_bin_secs);
    for i = 1:nbins 
        s1 = binned_trains{i}(1,:);
        s2 = binned_trains{i}(2,:);        
%         % smooth?  
%         win=hamming(5);t1smooth=conv(s1,win,'same');t2smooth=conv(s2,win,'same');
        pearson_xcov{i,1} = xcov(s1,s2,lag_units,'coef')';
        %LOW PASS 0.15%6th order, fc/fs/2 determined empiracally 
        %[b2,b1] = butter(6,0.03); pxcsmooth = filtfilt(b2,b1, pearson_xcov{i,1});
        %plot(pxcsmooth); hold on; plot(pearson_xcov{i,1},'y'); legend('smooth','none');
        %pearson_xcov{i,1} = pxcsmooth;
        count_xcor{i,1} = xcorr(s1,s2, lag_units)'; %how many times spike at same point
        %max(max([count_xcor{:,1}]));
        %time scale;
        time_scale=( (1:length(pearson_xcov{i,1})) - ((length(pearson_xcov{i,1})-1)/2) - 1)*params.time_bin_secs;
        pearson_xcov{i,2} = time_scale; count_xcor{i,3}= time_scale;
        %bin edges;
        pearson_xcov{i,3} = [bin_edges(i) bin_edges(i+1)];
        count_xcor{i,3}   = [bin_edges(i) bin_edges(i+1)];
        %plot(pearson_xcov); hold on%figure;plot(count_xcor)
    end
    pearson_xcov{1,3} = [-bin_edges(2) bin_edges(2)]; %correct for 0
    count_xcor{1,3}   = [-bin_edges(2) bin_edges(2)];
          
    fprintf('%.2f ',round(max(max(abs([pearson_xcov{:,1}]))), 2));
    
    
%GILAD STYLE
    %     smooth = 10e-5;%10e-3 %10e-5; %Empiracally so bins don't change rapidly for 8 bins!
%     tx = csaps( 1:length(c1.pt),double(c1.px), smooth);%10e-3; %10e-5;
%     ty = csaps( 1:length(c1.pt),double(c1.py), smooth);
%     dt = median(diff(c1.pt));
%     vx = fnval( fnder(tx),1:length(c1.pt))/dt; 
%     vy = fnval( fnder(ty),1:length(c1.pt))/dt;
    
        %spk1_by_bin = spk1_by_bin(~isnan(spk1_by_bin));  %remove nans
    %spk2_by_bin = spk2_by_bin(~isnan(spk2_by_bin)); %CHANGE INDICES????
    %assert(sum(isnan([ ... %{pos_by_bin;}% 
        %spk1_by_bin; spk2_by_bin]))==0,'nan bins compareByMovingDirection()');%+sum(isnan(spk_bins))
%     dir_bins_spk1 = cell(nbins,1);
%     dir_bins_spk2 = cell(nbins,1); 
%     for i = 1:length(spk1_by_bin)
%         if ~isnan(spk1_by_bin(i))
%             dir_bins_spk1{spk1_by_bin(i)}(end+1) = c1.st(i);
%         end
%     end
%     for i = 1:length(spk2_by_bin)
%         if ~isnan(spk2_by_bin(i))
%             dir_bins_spk2{spk2_by_bin(i)}(end+1) = c2.st(i);
%         end
%     end 
    


        %{
        t = cumsum(trains{i})
        S2=filter(ones(1,n)',1,M);
        S2=S2(n:n:end,:);
        padarray(mod something)
    
    
    figure(3); hold off; plot(0,0);%tx=c1.px(~isnan(pos_by_bin));ty=c1.py(~isnan(pos_by_bin));t=pos_by_bin(~isnan(pos_by_bin));
    for i = 1000:2000%length(pos_by_bin)
        if ~isnan(pos_by_bin(i)) && ~isnan(pos_by_bin(i-1))
           
        subplot(2,2,pos_by_bin(i)); hold on; title((pos_by_bin(i)-1)*90);
        plot(c1.px((i-1):i)-c1.px(i-1),c1.py((i-1):i)-c1.py(i-1),'b'); hold on
        %plot([0 vx(i)],[0 vy(i)],'r');
        end
        %plot(c1.px(i),c1.py(i),'x');
    end
    %}
%     gca.XAxisLocation = 'origin';
%     gca.YAxisLocation = 'origin';
%     
%     figure;
%     c = 10000;[vx(c) vy(c)]
%     subplot(2,2,1);quiver(0,0, vx(c),vy(c)); 
%     xlim([-20 20]);ylim([-20, 20])
%     subplot(2,2,2),plot([0 vx(c)],[0 vy(c)]);
%     hold on, plot(0,0,'x');xlim([-20 20]);ylim([-20, 20])
%     gca.XAxisLocation = 'origin';gca.YAxisLocation = 'origin';
%     subplot(2,2,3); hold on, plot(0,0,'x');
%     plot((c1.px( (c-1):c ) - c1.px(c-1)), (c1.py((c-1):c ) - c1.py(c-1)),'-');
%     xlim([-1 1]);ylim([-1, 1]);
%     subplot(2,2,4); hold off; plot(0,0,'x'); hold on,
%     quiver(20,0,vox(c),voy(c));
%     plot((c1.px( (c-1):c ) - c1.px(c-1))*40, 40*(c1.py((c-1):c ) - c1.py(c-1)),'-');    
%     xlim([-50 50]);ylim([-50, 50]);
%     axis('tight')
%     gca.XAxisLocation = 'origin';gca.YAxisLocation = 'origin';
    
    %{
    %smooth = 10e-5;
    %tx = fnval(csaps( 1:length(c1.pt),double(c1.px), smooth),1:length(c1.pt));
    %ty = fnval(csaps( 1:length(c1.pt),double(c1.py), smooth),1:length(c1.pt));
    %[b,a] = butter(6,0.1); %LOW PASS 0.03 0.15%6th order, fc/fs/2 determined empiracally
    %spx = filtfilt(b,a,double(c1.px));spy = filtfilt(b,a,double(c1.py));
    %%%
    pt = c1.pt; ox = c1.px; oy = c1.py; dtt=[0.02; diff(pt)];
    %svox = zeros(length(ox),1); vox = svox;
    %svoy = zeros(length(ox),1); voy = svoy;
    for i = 2:length(ox);
        %svox(i)= (spx(i)-spx(i-1))/(pt(i)-pt(i-1));svoy(i)= (spy(i)-spy(i-1))/(pt(i)-pt(i-1));
        vox(i)= (c1.px(i)-c1.px(i-1))/(pt(i)-pt(i-1));
        voy(i)= (c1.py(i)-c1.py(i-1))/(pt(i)-pt(i-1));
    end
    vox(1) = vox(2);voy(1) = voy(2);
    [b,a] = butter(6,0.06); %0.1 close
    svx = filtfilt(b,a,vox);svy = filtfilt(b,a,voy);
    spx = zeros(length(ox),1); spy = zeros(length(ox),1); spx(1) = ox(1); spy(1) = oy(1);
    for i = 2:length(ox);
        spx(i)= spx(i-1)+svx(i)*dtt(i);
        spy(i)= spy(i-1)+svy(i)*dtt(i);
    end
    %figure();plot(px);hold on; plot(vx);
    close all;
    for i = 1000:1000:10000
    s = i;d = 100;figure('Position', [0, 0, 2000, 1000]);  
    subplot(1,2,1); hold on; %plot(tx(s:s+d), ty(s:s+d),'r.'); plot(tx(s), ty(s),'rx'); 
    plot( ox(s:s+d),  oy(s:s+d),'go-');plot( ox(s), oy(s),'kx');
    plot(spx(s:s+d), spy(s:s+d),'ro-');plot(spx(s),spy(s),'kx');
    %plot(spx(s:s+d), spy(s:s+d),'mo-');plot(spx(s), spy(s),'mx');
    axis(gca,'tight'); title('x,y');
    subplot(1,2,2); hold on; %plot(vx(s:s+d), vy(s:s+d),'r.'); plot(vx(s), vy(s),'rx');
    %plot(vox(s:s+d), voy(s:s+d),'go-');plot(vox(s), voy(s),'gx');
    plot(svx(s:s+d), svy(s:s+d),'ro-');plot(svx(s), svy(s),'kx');axis(gca,'tight'); title('v');
    end
    subplot(3,2,3); plot(ox(s:s+d),'b.');title('x'); axis(gca,'tight'); set(gca,'xticklabel',[]); 
    subplot(3,2,4); plot(svx(s:s+d),'r.');title('vx'); axis(gca,'tight'); set(gca,'xticklabel',[]);
    subplot(3,2,5); plot(oy(s:s+d),'b.');title('y'); axis(gca,'tight'); set(gca,'xticklabel',[]);
    subplot(3,2,6); plot(svy(s:s+d),'r.');title('vy'); axis(gca,'tight'); set(gca,'xticklabel',[]);
    end
    %}
    %check offsets in position
    %figure;a = 0*10000+1; b = a+1000; plot(ox(a:b),oy(a:b)); hold on; plot(spx(a:b),spy(a:b));
    %plot([0 0], ylim,'black'),plot(xlim, [0 0],'black') 
    %plot([0 0], ylim,'black'),plot(xlim, [0 0],'black')

    %= sum(reshape(z(:,bin_size, ceil(length(z)/bin_size)));
    %z = sum(reshape([z;zeros(bin_size-mod(length(z),bin_size),1)],...
    %bin_size, ceil(length(z)/bin_size)));
    %binned_trains{i}(2,:) = z;
end

