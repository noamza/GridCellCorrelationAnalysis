function pcor0 = shuffleSpaceCorrelations (c1, c2, p)
    %[t1,t2] = createMsSpikeTrain(c1.st,100,c2.st); %WRONG because of aligning with xy in create ratemap
    t1=createMsSpikeTrain(c1.st); nb = p.nb;
    %offsets to deal with discrepency between 0 and first px py
    %so that during shifting we don't put spikes before px
    t2off=createMsSpikeTrain(c2.st-c2.pt(1)); %normalizes spike time relative to pt(1) = 1.
    t1smooth = movmean(t1,p.movmean);
    rm1 = createSmoothRateMap(c1,nb,t1smooth);
    %vars to speed up preallocation
    mx=max(c2.px); my=max(c2.py);
    rmt2 = histcounts2(c2.px,c2.py,0:mx/nb:mx,0:my/nb:my)';%save time
    t2smoff = zeros(1,ceil(c2.pt(1)*1000)+len(t2off)); %zeros array, length trial offset plus time of last spike. 
    st2 = (1:len(t2smoff))/1000;
    %shiftedi = zeros(1:len(t2off));
    %size(rm1)
    %delta to shift train by
    d = floor(len(t2off)/p.n);%right?  %amount to shift by each step
    pcor0 = []; 
    for i = 1:p.n %n+1 because wi = (i-1)*d; to include 0 shift
        %i
        %sprintf('done %.2f\%',round(i/(p.n+1),2));
        wi = (i-1)*d; %first one calculated is 0
        shiftedi = [wi+1:len(t2off) 1:wi]; %cyclic
        t2shifted = t2off(shiftedi); %shift offset spike train around by len/n
        t2smooth = movmean(t2shifted,p.movmean);
        %t2smooth = movmean(t2off([wi+1:len(t2off) 1:wi]),p.movmean);
        %t2smoff = [zeros(1,ceil(c2.pt(1)*1000)) t2smooth];
        t2smoff(ceil(c2.pt(1)*1000)+1:end) = t2smooth; %all zeros from 0 until pt(1) this is ok, 
        rm2 = createSmoothRateMap(c2,nb,t2smoff,rmt2,st2);%so times are always in range
        %scc = xcorr2(rm2-mean(rm2(:)),rm1-mean(rm1(:))); 
%         scc = normxcorr2(rm2,rm1); %reverse order for perspective should have negs
%         cenr = round(size(scc,2)/2); cenc = round(size(scc,1)/2);
%         scc0 = scc(cenr,cenc); %only check corr at 0,0 ??
%         pcor0(i) = scc0;
       scc0 = corrcoef(rm2,rm1); %same as norm2 of 0,0 but much faster
       pcor0(i) = scc0(2);
        %scc0
        %if mod(i-1,10) == 0
            %figure;imagesc(rm1); axis(gca,'square');hold on; plot(cenr, cenc,'mo');title(scc0); colormap jet;set(gca,'YDir','normal'); axis(gca,'tight'); set(gca,'Color','k');
            %figure;imagesc(rm2); axis(gca,'square');hold on; plot(cenr, cenc,'mo');title(scc0); colormap jet;set(gca,'YDir','normal'); axis(gca,'tight'); set(gca,'Color','k');
        %end
    end
end


%{ 

 tic
 i = 3; p.movmean = 25; p.n = 100;
 c1 = cells{pairs(i,1)}.before;
 c2 = cells{pairs(i,2)}.before;
 shuffleSpaceCorrelations (c1, c2, p);
 toc

 [m i] = max(abs(pcor0))

 q1 = c1.midall; q2 = c2.midall; 
 t1=createMsSpikeTrain(q1.st);
 t2=createMsSpikeTrain(q2.st);
 rm1 = createSmoothRateMap(q1,50,movmean(t1,25));
 rm2 = createSmoothRateMap(q2,50,movmean(t2,25));
 rm1 = q1.rm;
 rm2 = q2.rm;
 scc = xcorr2(rm2-mean(rm2(:)),rm1-mean(rm1(:))); %reverse order for perspective
 scc(round(size(scc,2)/2),round(size(scc,1)/2))
 figure;imagesc(rm1);figure;imagesc(rm2);figure;imagesc(scc);
%}