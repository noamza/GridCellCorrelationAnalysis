% MEAN FIRING RATE
fb=[];fd=[];fa=[];
for j=1:len(pairs)
    c1 = cells{pairs(j,1)}; c2 = cells{pairs(j,2)};
    s1 = c1.before.st; s1= len(s1)/(s1(end)-s1(1)); 
    s2 = c2.before.st; s2= len(s2)/(s2(end)-s2(1)); 
    fb(j,:) = [min(s1,s2) max(s1,s2) mean([s1 s2])];
    s1 = c1.midall.st; s1= len(s1)/(s1(end)-s1(1)); 
    s2 = c2.midall.st; s2= len(s2)/(s2(end)-s2(1));     
    fd(j,:) = [min(s1,s2) max(s1,s2) mean([s1 s2])];
    s1 = c1.after.st; s1= len(s1)/(s1(end)-s1(1)); 
    s2 = c2.after.st; s2= len(s2)/(s2(end)-s2(1));
    fa(j,:) = [min(s1,s2) max(s1,s2) mean([s1 s2])];
end
fa(fa(:,1)==inf,:)=[];

rsp=[];
for i=1:length(pairs)
    rsp(i,:) = [cells{pairs(i,1)}.midall.rayleigh_score,...
                cells{pairs(i,2)}.midall.rayleigh_score];
end
clhd = rsp(:,1)>0.4 & rsp(:,2)>0.4;
%FIRING RATE FIG

figure(1007); clf; xl = 'mean rate (Hz)'; yl = 'corr';
%BEFORE
% g score
subplot(231);d1 = fb(:,3); d2 = mb(:,3); plot(d1,d2,'.'); axis('tight');
title('PRE:: grid score');xlabel(xl); ylabel('score'); hold on;
f1 = fit(d1,d2,'poly1');plot(d1,f1(d1),'k-');
text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');
[];plot(d1(clhd), d2(clhd),'r.'); f1 = fit(d1(clhd), d2(clhd),'poly1'); plot(d1,f1(d1),'r-');
text(0.1,0.8,sprintf('hd: a=%.3f r=%.3f',round(f1.p1,3),round(ccof(d1(clhd),d2(clhd)),3)),...
    'Units','normalized','color','r');
f1 = fit(d1(~clhd), d2(~clhd),'poly1'); plot(d1,f1(d1),'b-');
text(0.1,0.7,sprintf('~hd: a=%.3f r=%.3f',round(f1.p1,3),round(ccof(d1(~clhd),d2(~clhd)),3)),...
    'Units','normalized','color','b');
%time
subplot(232);d1 = fb(:,3); d2 = ctsbma(:,1); plot(d1,d2,'.'); axis('tight'); 
title('temporal');xlabel(xl); ylabel(yl); hold on;
f1 = fit(d1,d2,'poly1');plot(d1,f1(d1),'k-');
text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');
[];plot(d1(clhd), d2(clhd),'r.'); f1 = fit(d1(clhd), d2(clhd),'poly1'); plot(d1,f1(d1),'r-');
text(0.1,0.8,sprintf('hd: a=%.3f r=%.3f',round(f1.p1,3),round(ccof(d1(clhd),d2(clhd)),3)),...
    'Units','normalized','color','r');
f1 = fit(d1(~clhd), d2(~clhd),'poly1'); plot(d1,f1(d1),'b-');
text(0.1,0.7,sprintf('~hd: a=%.3f r=%.3f',round(f1.p1,3),round(ccof(d1(~clhd),d2(~clhd)),3)),...
    'Units','normalized','color','b');
%space
subplot(233);d1 = fb(:,3); d2 = ctsbma(:,4); plot(d1,d2,'.'); axis('tight'); 
title('spatial');xlabel(xl); ylabel(yl); hold on;
f1 = fit(d1,d2,'poly1');plot(d1,f1(d1),'k-');
text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');
[];plot(d1(clhd), d2(clhd),'r.'); f1 = fit(d1(clhd), d2(clhd),'poly1'); plot(d1,f1(d1),'r-');
text(0.1,0.8,sprintf('hd: a=%.3f r=%.3f',round(f1.p1,3),round(ccof(d1(clhd),d2(clhd)),3)),...
    'Units','normalized','color','r');
f1 = fit(d1(~clhd), d2(~clhd),'poly1'); plot(d1,f1(d1),'b-');
text(0.1,0.7,sprintf('~hd: a=%.3f r=%.3f',round(f1.p1,3),round(ccof(d1(~clhd),d2(~clhd)),3)),...
    'Units','normalized','color','b');
%DURING
%score
subplot(234);d1 = fd(:,3); d2 = md(:,3); plot(d1,d2,'.'); axis('tight'); 
title('DUR:: grid score');xlabel(xl); ylabel('score'); hold on;
f1 = fit(d1,d2,'poly1');plot(d1,f1(d1),'k-');
text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');
[];plot(d1(clhd), d2(clhd),'r.'); f1 = fit(d1(clhd), d2(clhd),'poly1'); plot(d1,f1(d1),'r-');
text(0.1,0.8,sprintf('hd: a=%.3f r=%.3f',round(f1.p1,3),round(ccof(d1(clhd),d2(clhd)),3)),...
    'Units','normalized','color','r');
f1 = fit(d1(~clhd), d2(~clhd),'poly1'); plot(d1,f1(d1),'b-');
text(0.1,0.7,sprintf('~hd: a=%.3f r=%.3f',round(f1.p1,3),round(ccof(d1(~clhd),d2(~clhd)),3)),...
    'Units','normalized','color','b');
%time
subplot(235);d1 = fd(:,3); d2 = ctsbma(:,2); plot(d1,d2,'.'); axis('tight'); 
title('temporal');xlabel(xl); ylabel(yl); hold on;
f1 = fit(d1,d2,'poly1');plot(d1,f1(d1),'k-');
text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');
[];plot(d1(clhd), d2(clhd),'r.'); f1 = fit(d1(clhd), d2(clhd),'poly1'); plot(d1,f1(d1),'r-');
text(0.1,0.8,sprintf('hd: a=%.3f r=%.3f',round(f1.p1,3),round(ccof(d1(clhd),d2(clhd)),3)),...
    'Units','normalized','color','r');
f1 = fit(d1(~clhd), d2(~clhd),'poly1'); plot(d1,f1(d1),'b-');
text(0.1,0.7,sprintf('~hd: a=%.3f r=%.3f',round(f1.p1,3),round(ccof(d1(~clhd),d2(~clhd)),3)),...
    'Units','normalized','color','b');
%space
subplot(236);d1 = fd(:,3); d2 = ctsbma(:,5); plot(d1,d2,'.'); axis('tight'); 
title('spatial');xlabel('mean rate'); ylabel('corr'); hold on;
f1 = fit(d1,d2,'poly1');plot(d1,f1(d1),'k-');
text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');
[];plot(d1(clhd), d2(clhd),'r.'); f1 = fit(d1(clhd), d2(clhd),'poly1'); plot(d1,f1(d1),'r-');
text(0.1,0.8,sprintf('hd: a=%.3f r=%.3f',round(f1.p1,3),round(ccof(d1(clhd),d2(clhd)),3)),...
    'Units','normalized','color','r');
f1 = fit(d1(~clhd), d2(~clhd),'poly1'); plot(d1,f1(d1),'b-');
text(0.1,0.7,sprintf('~hd: a=%.3f r=%.3f',round(f1.p1,3),round(ccof(d1(~clhd),d2(~clhd)),3)),...
    'Units','normalized','color','b');

suptitle('Mean Firing Rate (row 1 before, row 2 during)');

a = [];
for j=1:len(cels)
    c = cells{j}.before;
    a(j,:) = [1/(min(diff(c.st))),c.max_r]; %max rate
end

for j=1:len(cels)
c = cells{cels(j)}.before; b=1./diff(c.st);
[sum(b<1) sum(b>1) sum(b>100) sum(b>500) mean(b(b<500 & b>1))]
end
