function pcor0 = spaceCorrelationSmoothed(c1,c2,p)

    t1=createMsSpikeTrain(c1.st); nb = p.nb;
    t1s = movmean(t1,p.movmean);
    rm1 = createSmoothRateMap(c1,nb,t1s);

    t2=createMsSpikeTrain(c2.st);
    t2s= movmean(t2,p.movmean);
    rm2 = createSmoothRateMap(c2,nb,t2s);

    %     scc = normxcorr2(rm2,rm1); %reverse order for perspective should have negs
    %     sccc = corrcoef(rm2,rm1);
    %     cenr = round(size(scc,2)/2); cenc = round(size(scc,1)/2);
    %     scc0 = scc(cenr,cenc); %only check corr at 0,0 ??
    %     pcor0 = [scc0 sccc(2)];%corr(rm2-rm1)];
    pcor0=ccof(rm2,rm1);
    
%     figure()
%     subplot(211); imagesc(rm1);
%     subplot(212); imagesc(rm2);

end