%SUPP SCORE VS CORR 
mb=[];md=[];
for j=1:len(pairs)
    c1 = cellsn(pairs(j,1)); c2 = cellsn(pairs(j,2));%c1.ind, c2.ind
    s1 = c1.before.gridscore; s2 = c2.before.gridscore; s1(s1==-2)=0; s2(s2==-2)=0;
    mb(j,:) = [min(s1,s2) max(s1,s2) mean([s1 s2])];
    s1 = c1.midall.gridscore; s2 = c2.midall.gridscore; s1(s1==-2)=0; s2(s2==-2)=0;
    md(j,:) = [min(s1,s2) max(s1,s2) mean([s1 s2])];
end


figure(1006); clf; tos=1;

subplot(241);d1 = mb(:,3); d2 = ctsbma(:,tos); plot(d1,d2,'.'); axis('tight'); 
title('time corr pre');xlabel('score'); ylabel('corr'); hold on;
f1 = fit(d1,d2,'poly1');plot(d1,f1(d1),'-');
text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');

subplot(242);d1 = md(:,3); d2 = ctsbma(:,tos+1); plot(d1,d2,'.'); axis('tight'); 
title('time corr dur');xlabel('score'); ylabel('corr'); hold on;
f1 = fit(d1,d2,'poly1');plot(d1,f1(d1),'-');
text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');

tos = 4;

subplot(243);d1 = mb(:,3); d2 = ctsbma(:,tos); plot(d1,d2,'.'); axis('tight'); 
title('space corr pre');xlabel('score'); ylabel('corr'); hold on;
f1 = fit(d1,d2,'poly1');plot(d1,f1(d1),'-');
text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');

subplot(244);d1 = md(:,3); d2 = ctsbma(:,tos+1); plot(d1,d2,'.'); axis('tight'); 
title('space corr dur');xlabel('score'); ylabel('corr'); hold on;
f1 = fit(d1,d2,'poly1');plot(d1,f1(d1),'-');
text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');

tos=1;

subplot(245);d1 = mb(:,3); d2 = ctsbma(:,tos); d2 = abs(d2); plot(d1,d2,'.'); axis('tight'); 
title('time corr pre');xlabel('score'); ylabel('abs(corr)'); hold on;
f1 = fit(d1,d2,'poly1');plot(d1,f1(d1),'-');
text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');

subplot(246);d1 = md(:,3); d2 = ctsbma(:,tos+1); d2 = abs(d2); plot(d1,d2,'.'); axis('tight'); 
title('time corr dur');xlabel('score'); ylabel('abs(corr)'); hold on;
f1 = fit(d1,d2,'poly1');plot(d1,f1(d1),'-');
text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');

tos = 4;

subplot(247);d1 = mb(:,3); d2 = ctsbma(:,tos); d2 = abs(d2); plot(d1,d2,'.'); axis('tight'); 
title('space core pre');xlabel('score'); ylabel('abs(corr)'); hold on;
f1 = fit(d1,d2,'poly1');plot(d1,f1(d1),'-');
text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');

subplot(248);d1 = md(:,3); d2 = ctsbma(:,tos+1); d2 = abs(d2); plot(d1,d2,'.'); axis('tight'); 
title('space core dur');xlabel('score'); ylabel('abs(corr)'); hold on;
f1 = fit(d1,d2,'poly1');plot(d1,f1(d1),'-');
text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');

suptitle('Mean Grid Score vs Correlation');