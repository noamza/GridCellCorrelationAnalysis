% MEAN FIRING RATE
%load('C:\Noam\Data\muscimol\cells15nan')
load('.\data\shuffling1000nanv2')
load('.\data\pairsv2')
fb=[];fd=[];fa=[]; rst=0.4; %NEED TO VE CONSISTENT
gridscore='gridscore'; 
%for v3%
%gridscore='gs2';

for j=1:len(pairs)
    c1 = cellsn(pairs(j,1)); c2 = cellsn(pairs(j,2));
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

gb=[];gd=[];
for j=1:len(pairs)
    c1 = cellsn(pairs(j,1)); c2 = cellsn(pairs(j,2));%c1.ind, c2.ind
    s1 = c1.before.(gridscore); s2 = c2.before.(gridscore); s1(s1==-2)=0; s2(s2==-2)=0;
    gb(j,:) = [min(s1,s2) max(s1,s2) mean([s1 s2])];
    s1 = c1.midall.(gridscore); s2 = c2.midall.(gridscore); s1(s1==-2)=0; s2(s2==-2)=0;
    gd(j,:) = [min(s1,s2) max(s1,s2) mean([s1 s2])];
end

rsp=[];
for i=1:length(pairs)
    rsp(i,:) = [cellsn(pairs(i,1)).midall.rayleigh_score,...
                cellsn(pairs(i,2)).midall.rayleigh_score];
end
clhd = rsp(:,1)>rst & rsp(:,2)>rst; 

fig=figure(1007); clf; set(fig,'color','w', 'Position', [200 70 1100 550]);
arg=[];arg.setlim=false; yl='corr';
%BEFORE
% g score
subplot(231);x = fb(:,3); y = gb(:,3); fp(x,y,clhd); 
title('pre grid score');ylabel('score'); 
%time
subplot(232);x = fb(:,3); y = ctsbma(:,1); fp(x,y,clhd); 
title('pre temporal');ylabel(yl); 

%space
subplot(233);x = fb(:,3); y = ctsbma(:,4);fp(x,y,clhd); 
title('pre spatial');ylabel(yl);
%DURING
%score
subplot(234);x = fd(:,3); y = gd(:,3); fp(x,y,clhd); 
title('dur grid score');ylabel('score');

%time
subplot(235);x = fd(:,3); y = ctsbma(:,2);fp(x,y,clhd); 
title('dur temporal');ylabel(yl); 

%space
subplot(236);x = fd(:,3); y = ctsbma(:,5);fp(x,y,clhd); 
title('dur spatial');ylabel('corr'); 

suptitle('mean firing rate vs grid score and correlation');

function fp(x,y,clhd)
plot(x,        y       ,'k.'); hold on;
plot(x( clhd), y( clhd),'r.');
plot(x(~clhd), y(~clhd),'b.');
f1 = fit(x,        y,'poly1');        plot(x,f1(x),'k-');
f2 = fit(x(clhd),  y( clhd),'poly1'); plot(x,f2(x),'r-');
f3 = fit(x(~clhd), y(~clhd),'poly1'); plot(x,f3(x),'b-');
axis tight; axis square;
text(0.1,0.9,sprintf('all: a=%.3f r=%.3f',round(f1.p1,3),round(ccof(x,       y)       ,3)),...
    'Units','normalized');
text(0.1,0.8,sprintf('hd : a=%.3f r=%.3f',round(f2.p1,3),round(ccof(x(clhd), y(clhd)) ,3)),...
    'Units','normalized','color','r');
text(0.1,0.7,sprintf('!hd: a=%.3f r=%.3f',round(f3.p1,3),round(ccof(x(~clhd),y(~clhd)),3)),...
    'Units','normalized','color','b');
xlabel('firing rate (Hz)');
end

%{
a = [];
for j=1:len(cels)
    c = cellsn(j).before;
    a(j,:) = [1/(min(diff(c.st))),c.max_r]; %max rate
end

for j=1:len(cels)
c = cellsn(cels(j)).before; b=1./diff(c.st);
[sum(b<1) sum(b>1) sum(b>100) sum(b>500) mean(b(b<500 & b>1))]
end
%}
