function gshuf = shuffleGridscoreNan(c,nshuf,nb,movm,sig,pval,s)  
    rt=histcounts2(c.px,c.py, 0:max(c.px)/nb:max(c.px),0:max(c.py)/nb:max(c.py))';
    sta = (1:ceil(c.pt(end)*1000))/1000;
    d = (c.pt(end)-c.pt(1))/(nshuf); %how much to shift by %/(nshuf+1) 
    gshuf = zeros(nshuf,1);
    if isequal(s,'midall') %for checking failure
        %gshuf = zeros(nshuf,1)+3;
    end
    i=0;
    while i < nshuf % <= %for /(nshuf+1) 
        i=i+1;
        wi = (i-1)*d; %first one calculated is 0  
        t = c.st + wi; t(t>c.pt(end)) = t(t>c.pt(end))-c.pt(end)+c.pt(1); %shift cells to beginning
        tms = createMsSpikeTrain(t);
        tsm = movmean(tms,movm);
        if length(tsm) <= length(sta)
            st=sta(1:length(tsm)); %should be:  st = (1:len(trs))/1000;
        else
            st=(1:length(tsm))/1000;
        end
        rm = createSmoothRateMapNan(c,nb,tsm,rt,st);
        ac = xcorr2g(rm,rm);
        %figure(1);imgsc(ac);title(i);
        gshuf(i)=gridscore2(ac,sig);
        
%         if i==1
%             '[gshuf(i)  c.gridscore]'
%              [gshuf(i)  c.gridscore]
%         end
        
        if isequal(s,'midall') % fail if too significant
             % positive gridscore is greater than p percent of shuffled scores % IS significant
            if gshuf(1) >= 0.7 || gshuf(1) <= 0 ||...
               sum( gshuf(2:end) <= gshuf(1) ) >= ceil(nshuf*pval) %&& gshuf(1) > 0
                ['mid:TOO significant gridscore shuffling pval=' n2(pval)]
                ['gshuf(1) ' n2(gshuf(1))]
                i=inf;
            end
        else % before or after, fail if not significant 
             % number of shuffled gridscores is greater than p percent %NOT significant
            if sum( gshuf(2:end) > gshuf(1) ) >= ceil(nshuf*pval)  %|| gshuf(1) <= 0
                ['bef:NOT significant gridscore shuffling i=' n2(i) ' pval=' n2(pval) ' more than ' n2(ceil(nshuf*pval))]
                ['gshuf(1) ' n2(gshuf(1))]
                i=inf;
            end
        end
    end
end



