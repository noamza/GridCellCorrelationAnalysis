function pcor0 = shuffleSpace2Correlations (c1, c2, p)    
    t1=createMsSpikeTrain(c1.st);
    
    t1s = movmean(t1,p.movmean);
    
    rm1 = createSmoothRateMapNan(c1,p.nb,t1s);%createSmoothRate()
    
    rt2=histcounts2(c2.px,c2.py, 0:max(c2.px)/p.nb:max(c2.px),0:max(c2.py)/p.nb:max(c2.py))';
    
    st = (1:ceil(c2.pt(end)*1000))/1000;
    
    d = (c2.pt(end)-c2.pt(1))/p.n;%right?  %amount to shift by each step
    %b=a+d; b(b>c(end))=b(b>c(end))-c(end)+c(1)    
      
    pcor0 = []; 
    for i = 1:p.n
        
        wi = (i-1)*d; %first one calculated is 0  
        t = c2.st + wi; t(t>c2.pt(end)) = t(t>c2.pt(end))-c2.pt(end)+c2.pt(1); %shift cells to beginning
        t2 = createMsSpikeTrain(t);
        %[t1,t2] = removeOverlappingSpikes(t1,t2, p.remove); %USE??
        % t1=createMsSpikeTrain(t1);       
        t2s = movmean(t2,p.movmean);
        st2=st(1:len(t2s)); %should be:  st = (1:len(trs))/1000;
        rm2 = createSmoothRateMapNan(c2,p.nb,t2s,rt2,st2);  %createSmoothRate()
        pcor0(i)=ccof(rm2,rm1);
    end
end

%         if i == 1 || i == 50
%             'here'
%             figure
%             imagesc(rm2);
%         end