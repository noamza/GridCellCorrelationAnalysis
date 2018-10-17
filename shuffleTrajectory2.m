function [ output_args ] = shuffleTrajectory2( input_args )


pairs = [2,7;11,12;11,13;11,14;12,13;12,14;13,14;22,24;22,28;24,28;34,35;34,37;34,39;34,41;35,37;35,39;35,41;37,39;37,41;39,41;51,53;61,62;61,66;61,68;62,66;62,68;66,68;77,78;77,79;78,79;80,81;82,86;100,101;100,103;100,104;100,105;100,111;101,103;101,104;101,105;101,111;103,104;103,105;103,111;104,105;104,111;105,111;112,114;112,115;112,116;112,117;112,118;112,119;114,115;114,116;114,117;114,118;114,119;115,116;115,117;115,118;115,119;116,117;116,118;116,119;117,118;117,119;118,119;171,173;171,174;171,176;171,177;173,174;173,176;173,177;174,176;174,177;176,177;178,179;178,180;178,182;178,184;179,180;179,182;179,184;180,182;180,184;182,184;188,189;199,201;199,209;201,209;230,231;274,276;274,278;276,278;282,283;282,287;283,287];

b = []; m = []; a = [];

h = figure(3534);c = cells{1}.before; tr = createMsSpikeTrain(c.st);
nbins = 300; nmap = jet; nmap(1,:) = [0 0 0];
for i = 1:20%length(pairs{    
    %c = cells{pairs(i,1)}; c = c.before;
    %tr = makeMsSpikeTrain(c.st,100);
    mx = max(c.px);
    my = max(c.py);
    pxi = discretize(c.px, 0:mx/nbins:mx);
    pyi = discretize(c.py, 0:my/nbins:my);
    %t = diff(c.pt(pw)); t = [median(t); t];
    rmt = accumarray([pyi pxi], 1, [nbins nbins]);
    %rmt = accumarray([pyi pxi], t, [nbins nbins]);
    %rmt(rmt<1) = 1e10;
    win = 2^(i-1);
    trs = movmean(tr,win);
    st = [1:len(trs)]/1000;
    sti = discretize(st, [-Inf, mean([c.pt(2:end)'; c.pt(1:end-1)']), +Inf]);
    sxi = discretize(c.px(sti), 0:mx/nbins:mx);
    syi = discretize(c.py(sti), 0:my/nbins:my);
    rmss = accumarray([syi sxi], trs, [nbins nbins]);
    rm2 = rmss./rmt;
    rm2(isnan(rm2)) = 0;
    %subplot(5,4,i);
    ax= axes('position',apos(4,20,i));
    imagesc(imgaussfilt(rm2,1));colormap(nmap); axis(gca,'square');
    %plot(trs)
    set(gca,'YDir','normal'); axis(gca,'tight'); set(gca,'Color','k');
    title(sprintf('win %.dms',win),'color','m','fontsize', 12);
end
%HD HD HD HD
h = figure(5534);c = cells{6}.before; tr2 = createMsSpikeTrain(c.st); rbin = linspace(-pi,pi,1000);
sti = discretize([1:len(trs)]/1000, [-Inf,mean([c.pt(2:end)'; c.pt(1:end-1)']), +Inf]);
c.si = discretize(c.st, [-Inf, mean([c.pt(2:end)'; c.pt(1:end-1)']), +Inf]);
for i = 1:20%length(pairs{    
    win = 2^(i-1);
    trs = movmean(tr2,win);
    shdi = discretize(wrapToPi(c.hd(sti)), rbin);
    rd2 = accumarray(shdi, trs,[len(rbin) 1])';
    rd = accumarray(discretize(wrapToPi(c.hd(c.si)), rbin), 1,[len(rbin) 1])';
    rdt = accumarray(discretize(wrapToPi(c.hd), rbin), 1, [len(rbin) 1])';
    rd2a = rd2./rdt;  rda = rd./rdt;
    rd2a(isnan(rd2a)) = 0; rda(isnan(rda)) = 0;
    r2 = sqrt( sum(cos(rbin).*rd2a).^2+sum(sin(rbin).*rd2a).^2)/sum(rd2a);
    ax= axes('position',apos(4,20,i));
    plot(ax, rbin,rda);
    hold on;
    plot(ax, rbin,rd2a);  axis(ax,'tight');
    title(ax,sprintf('win %.dms rayleigh orig=%.2f smth=%.2f',win,c.rayleigh_score,r2),'fontsize',9);
end


     nbins = 50;
        mx = max(c.px);
        my = max(c.py);
        pxi = discretize(c.px, 0:mx/nbins:mx);
        pyi = discretize(c.py, 0:my/nbins:my);
        %t = diff(c.pt(pw)); t = [median(t); t];
        rmt = accumarray([pyi pxi], 1, [nbins nbins]);
        %rmt = accumarray([pyi pxi], t, [nbins nbins]);
        rmt(rmt<1) = 1e10;
        %rmt(rmt<min(t)) = 1e10;
        %rmt = accumarray([pxi' pyi'], 1, [nbins nbins]); %SORT OUT WITH FUNCTION
        sxi = discretize(c.px(sw), 0:mx/nbins:mx);
        syi = discretize(c.py(sw), 0:my/nbins:my);
        rmss = accumarray([syi sxi], 1, [nbins nbins]);
        %rms = accumarray([sxi' syi'], 1, [nbins nbins]);
        rm = rmss./(rmt*dt);
        %rm = rms./rmt;
        [m, i] = max(rm(:));
        ind2sub(size(rm),i);
        %t = rms-rmt; t(t<0)=0; rmt = rmt + t; %add in extra timestep ONLY when spikes occured faster than timestep
        rm = (rmss ./ (rmt*dt));
        rm(isnan(rm)) = 0; %DO THIS?
        %maxfiring rate
        %rm = imgaussfilt(rm,1);
        %rm = rm/max(rm(:)); %normalize %REMOVE??        
        %f = gaussian2d(nbins,2);
        %imagesc(ax, conv2(rm,f,'same'));
        figure(343);
        imagesc(ax, imgaussfilt(c.rm,1)); %rm
        title(ax,sprintf('max %.1fHz',m),'color','m','fontsize', tfs);
        
    end
end



end

