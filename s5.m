%EXTENDED TIME WINDOW
load('.\data\shuffling1000nanV2')
strt = 15 *60; midt = (45/1) *60; endt=45 *60 +inf;
fig=figure(1005); clf; set(fig,'color','w', 'Position', [200 70 700 500]);
xl = 'pre'; yl = 'dur';
%TIME
%1st half
subplot(221);x = ctsbma(:,1); y = wd(:,1);arg=[];
plotARP(x,y,arg);xlabel(xl); ylabel(yl);
title(sprintf('[15+ : %.0f min]',midt/60));
%2nd half
subplot(222);x = ctsbma(:,1); y = wd(:,2);plotARP(x,y,arg);
xlabel(xl); ylabel(yl);
title(sprintf('[%.0f : %.0f min]',midt/60,endt/60));

%SPACE
subplot(223);x = ctsbma(:,4); y = wd(:,3);
plotARP(x,y,arg);xlabel(xl); ylabel(yl);
%2nd half
subplot(224);x = ctsbma(:,4); y = wd(:,4);
plotARP(x,y,arg);xlabel(xl); ylabel(yl);

suptitle('correlations by muscimol time window (row 1 temporal, row 2 spatial)');


%{
%load('C:\\Noam\\Data\\muscimol\\noam\\cells_Infmin_d_patchtraj_rayleigh'); %personal

tic; wd=[]; wb=[]; %was 100 (anylysis done in 50?)
for j=1:len(pairs)
    j
    c1 = cells{pairs(j,1)}; c2 = cells{pairs(j,2)};
    ['midall end ' n2(c1.midall.pt(end)/60)]
    
    p=[]; p.movmean = 25; nbins = 50; %spatial
%     wb(j) = c1.midall.pt(end)-c1.midall.pt(1);
%    t = c1.midall.pt; tm = t(1)+(t(end)-t(1))/2;
%     [train1,train2] = createMsSpikeTrain(c1.midall.st,t(end),c2.midall.st, t(1),tm);
%     [a1, ~, t1sm, t2sm] = timeCorrelationSmoothed( train1,train2,p);
    
    %im = ceil(len(c1.midall.pt)/2); find index of half way
    strt = 15 *60; midt = (45/1) *60; endt=45 *60 +inf; assert nothing starts before 15
    %is = find(c1.midall.pt>=strt*60,1);
    %im = find(c1.midall.pt>=midt*60,1);
    %ie = len(c1.midall.pt);
    %ie = find(c1.midall.pt<=endt*60,1,'last');
    
    
    %c1a = windowsesh(c1.midall,is,im); c2a = windowsesh(c2.midall,is,im); %was 1
    %c1b = windowsesh(c1.midall,im,ie); c2b = windowsesh(c2.midall,im,ie);
    c1a = windowsesh(c1.midall,0,0,strt,midt); c2a = windowsesh(c2.midall,0,0,strt,midt); %was 1
    c1b = windowsesh(c1.midall,0,0,midt,endt); c2b = windowsesh(c2.midall,0,0,midt,endt);
    
    if ~(isempty(c1a) || isempty(c1b) || isempty(c2a) || isempty(c2b) )
        %1st half
        [t1,t2] = createMsSpikeTrain(c1a.st,-1,c2a.st);
        [ta, ~, t1sm, t2sm] = timeCorrelationSmoothed(t1,t2,p);
        sa=ccof(createSmoothRateMapNan(c2a,nbins,t2sm),createSmoothRateMap(c1a,nbins,t1sm));
        %2nd half
        [t1,t2] = createMsSpikeTrain(c1b.st,-1,c2b.st);
        [tb, ~, t1sm, t2sm] = timeCorrelationSmoothed(t1,t2,p);
        sb=ccof(createSmoothRateMapNan(c2b,nbins,t2sm),createSmoothRateMap(c1b,nbins,t1sm));
    else
        ta=0;tb=0;sa=0;sb=0;
    end
    sa(sa==1)=0;sb(sb==1)=0;
    wd(j,:) = [ta tb sa sb];  
end
wd(isnan(wd))=0; toc
%}

