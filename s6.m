%load('C:\Noam\Data\muscimol\cells15nan')
load('.\data\shuffling1000nanv2')
load('.\data\pairsv2')
gridscore='gridscore'; 
%for v3%
%gridscore='gs2';

%GRIDSCORE VS CORR
mb=[];md=[];
for j=1:len(pairs)
    c1 = cellsn(pairs(j,1)); c2 = cellsn(pairs(j,2));%c1.ind, c2.ind
    s1 = c1.before.(gridscore); s2 = c2.before.(gridscore); s1(s1==-2)=0; s2(s2==-2)=0;
    mb(j,:) = [min(s1,s2) max(s1,s2) mean([s1 s2])];
    s1 = c1.midall.(gridscore); s2 = c2.midall.(gridscore); s1(s1==-2)=0; s2(s2==-2)=0;
    md(j,:) = [min(s1,s2) max(s1,s2) mean([s1 s2])];
end


fig=figure(1006); clf; set(fig,'color','w', 'Position', [200 70 1100 550]);
arg=[];arg.setlim=false; xl='score'; yl='corr'; 

tos=1;%time pre
subplot(241);x = mb(:,3); y = ctsbma(:,tos); plotARP(x,y,arg);
title('time pre');xlabel(xl); ylabel(yl);
subplot(242);x = md(:,3); y = ctsbma(:,tos+1);plotARP(x,y,arg);
title('time dur');xlabel(xl); ylabel(yl); 

tos = 4;%spatial pre
subplot(243);x = mb(:,3); y = ctsbma(:,tos); plotARP(x,y,arg);
title('spatial pre');xlabel(xl); ylabel(yl); 
subplot(244);x = md(:,3); y = ctsbma(:,tos+1);plotARP(x,y,arg);
title('spatial dur');xlabel(xl); ylabel(yl); 

tos=1; yl = 'abs(corr)';%time pre
subplot(245);x = mb(:,3); y = ctsbma(:,tos); y = abs(y); plotARP(x,y,arg);
title('time pre');xlabel(xl); ylabel(yl); 
subplot(246);x = md(:,3); y = ctsbma(:,tos+1); y = abs(y); plotARP(x,y,arg); 
title('time dur');xlabel(xl); ylabel(yl); 

tos = 4; %spatial pre
subplot(247);x = mb(:,3); y = ctsbma(:,tos); y = abs(y); plotARP(x,y,arg);
title('spatial pre');xlabel(xl); ylabel(yl); 
subplot(248);x = md(:,3); y = ctsbma(:,tos+1); y = abs(y); plotARP(x,y,arg); 
title('spatial dur');xlabel(xl); ylabel(yl); 

suptitle('Mean Grid Score vs Correlation');