%{
% Load before calling function!!
load('.\data\shuffling1000nan');
load('.\data\pairs'); cels = unique(pairs(:))';
%}

function f4(cellsn, pairs)



%% preamble
dbstop if error 

load('.\\data\\shuffling1000nan','ctsbma','ptsbma','pptsbma');

fig = figure(994); 
set(fig,'color','w', 'Position', [200 70 1200 900]);

fss='fontsize';mss='markersize';p='Parent';fwss='fontweight';bss='bold';vss='visible';noss='none';
bcss='BackgroundColor'; btss='BorderType';spss='spacing';dmss='DividerMarkings';ofss='off';onss='on';
hzss='HorizontalAlignment';cnss='center'; coss='color';unss='Units';nzss='normalized';
s = {'before';'midall';'after'};tss = {'pre','dur','post'};
vs='rayleigh_score';va='rayleigh_angle';
%fs ={fss 14}; afs ={fss 10}; mfs ={mss 10}; tfs = {fss,12}; fw={fwss,bss};
tx=0;ty=0.5; ty2=0.9;
fs = 12; afs = 10; mfs = 10; tfs = 12;  lfs = 3; h = 25;

rst=0.45; rsth=rst; rstl=rst;
ci = [112, 114, 117, 118, 119]; %ci = [103, 104, 105, 100,111];  ci = [282,283,287]; ci = [112, 114, 115, 116, 117,118,119]; 115
cels = unique(pairs(:))';
clls = cellsn(cels); %t = [clls.midall]; t = [t.rayleigh_score]; hdclls = clls(t>rst); t=[];



%% make clusters
rsp=[];rs=[]; rmi=true(len(ctsbma),1);
for i=1:length(pairs)
    rsp(i,:) = [cellsn(pairs(i,1)).before.rayleigh_score,...
                cellsn(pairs(i,2)).before.rayleigh_score,...
                cellsn(pairs(i,1)).midall.rayleigh_score,...
                cellsn(pairs(i,2)).midall.rayleigh_score];
end
for i=1:length(cels)
    rs(i,:)  = [cellsn(cels(i)).before.rayleigh_score,...
        cellsn(cels(i)).midall.rayleigh_score];
end
ccl2 = rs(:,2)>=rsth & rs(:,1)<rstl; 
ccl1 = rs(:,2) <rstl & rs(:,1)<rstl; %ccl1 = ~ccl2; 
%rsp: rscore[c1b c2b c1m c2m]
cl2 = rsp(:,3)>=rsth & rsp(:,4)>=rsth ...   %midall
    & rsp(:,1)< rstl & rsp(:,2)< rstl;      %before
cl1 = rsp(:,3)< rstl & rsp(:,4)< rstl ...     %midall
    & rsp(:,1)< rstl & rsp(:,2)< rstl;     %before
%cl1 = ~cl2;
chd = cels(ccl2);
hdclls=cellsn(chd);
[sum(cl2) sum(cl1)];
rmi = ctsbma(:,1)<0.15 & ctsbma(:,2)<0.07;

%% gTop
gTop = uix.GridFlex(p, fig,spss,5, bcss,'w',dmss,ofss);
%% gLeft
gLeft = uix.GridFlex(p, gTop,spss,5, bcss,'w',dmss,ofss);

%% A
%axes(p,uicontainer(p,gleft,bcss,'w'),vss,ofss);
axes(p,uix.Panel(p,gLeft,bcss,'w',btss,noss),vss,ofss);
text(tx,ty,'A',fwss,bss,fss,fs); axis off;
gA = uix.GridFlex(p,gLeft,spss,1, bcss,'w',dmss,onss);
%%{
axA = {}; pf = []; s1 = [];
for i = 1:len(ci)
    c = cellsn(ci(i));
    %TRAJ
    axes(p,uicontainer(p,gA,bcss,'w'),vss,ofss);
    text(0.5,0.5,sprintf('C%d',i),fwss,bss,fss,mfs,hzss,cnss);
            s1.x=''; s1.y = s1.x; s1.t ='';
    axA{end+1} = axes(p,uicontainer(p,gA,bcss,'w'));
    axA{end}   = plotTR(axA{end},c.before,s1); s1.t = '';
    axA{end+1} = axes(p,uicontainer(p,gA,bcss,'w'));
    axA{end}   = plotTR(axA{end},c.midall,s1);
    %AC
            s1.t = sprintf('grid=%.1f',round(c.before.gridscore,2));
    axA{end+1} = axes(p,uicontainer(p,gA,bcss,'w'));
    axA{end}   = plotAC(axA{end},c.before,s1);
            s1.t = sprintf('grid=%.1f',round(c.midall.gridscore,2));
    axA{end+1} = axes(p,uicontainer(p,gA,bcss,'w'));
    axA{end}   = plotAC(axA{end},c.midall,s1);
    %RAY
            s1.t = sprintf('rayleigh=%.1f',round(c.before.rayleigh_score,2));
    axA{end+1} = axes(p,uicontainer(p,gA,bcss,'w'));
    axA{end}   = plotHD(axA{end},c.before,s1,1);
            s1.t = sprintf('rayleigh=%.1f',round(c.midall.rayleigh_score,2));
    axA{end+1} = axes(p,uicontainer(p,gA,bcss,'w'));
    axA{end}   = plotHD(axA{end},c.midall,s1,1);
end

set(gA,'Heights', [15 -1 -1 -1 -1 -1 -1]);%,'Widths', [-1 -1 -1]);
%}

%% B
axes(p,uix.Panel(p,gLeft,bcss,'w',btss,noss),vss,ofss);
text(0,0.0,'B',fwss,bss,fss,fs,hzss,cnss);  axis off;
gB = uix.GridFlex(p,gLeft,spss,5, bcss,'w',dmss,ofss);
%%{
axB = {}; pf = []; s1 = [];
s1.x=''; s1.y = ''; s1.t = ''; pf.off = 0; pf.lag = 1000; pf.movmean = 100;
for i = 1:len(ci)-1
    for ii = i+1:len(ci)
        c1 = cellsn(ci(i));c2 = cellsn(ci(ii));
        if i+ii>3
            %uix.Empty(p,gB);
        end
        %s1.t = sprintf('C%dxC%d',i,ii);
        axes(p,uicontainer(p,gB,bcss,'w'),vss,ofss);
        text(0.3,0.5,sprintf('C%dxC%d',i,ii),fwss,bss,fss,mfs,hzss,cnss);
        axB{end+1} = axes(p,uicontainer(p,gB,bcss,'w'));
        axB{end}  = plotTimeCorr(axB{end},c1.before,c2.before,pf,s1);
        axB{end+1} = axes(p,uicontainer(p,gB,bcss,'w'));
        s1.t = '';
        axB{end}  = plotTimeCorr(axB{end},c1.midall,c2.midall,pf,s1); 
        axB{end+1} = axes(p,uicontainer(p,gB,bcss,'w')); 
        axB{end}  = plotSpaceCorr(axB{end},c1.before,c2.before,pf,s1);
        axB{end+1} = axes(p,uicontainer(p,gB,bcss,'w')); 
        axB{end}  = plotSpaceCorr(axB{end},c1.midall,c2.midall,pf,s1);
        
    end
end
set(gB,'Heights', [15 -1 -1 -1 -1]);%,'Widths', [-1 -1 -1]);
%}

set(gLeft,'Heights', [25 -2 25 -1]);

%gCenter
gCent = uix.GridFlex(p, gTop,spss,5, bcss,'w',dmss,ofss);

%EMPTY
axes(p,uicontainer(p,gCent,bcss,'w'),vss,ofss);
plot(1,1,'mo'); hold on %,mss,mfs
plot(1,1,'b.',mss,mfs);
plot(1,1,'r.',mss,mfs);
plot(1,1,'w.',mss,mfs+15);box off; axis off;
legend({'cohort';'cluster 1';'cluster 2'},fss,afs,'Location','south');

%% C hd single
gC = uix.GridFlex(p,gCent,spss,1, bcss,'w',dmss,ofss);
axes(p,uicontainer(p,gC,bcss,'w'),vss,ofss);
text(0,0.9,'C',fwss,bss,fss,fs);
axes(p,uicontainer(p,gC,bcss,'w'),vss,ofss);
plot(rs(:,1),rs(:,2),'mo'); hold on;  %,mss,mfs
%hold on; plot(rs(rs(:,2)>0.5,1),rs(rs(:,2)>0.5,2),'r.',mss,mfs);
plot(rs(ccl1,1),rs(ccl1,2),'b.',mss,mfs);
plot(rs(ccl2,1),rs(ccl2,2),'r.',mss,mfs);
%title('rayleigh score by cell');
%legend({'cohort';'cluster 1';'cluster 2'},fss,afs,'Location','east');
set(gca,'XTick',[tx,ty,1]);set(gca,'XTickLabel',[tx,ty,1],fss,afs)
set(gca,'YTick',[tx,ty,1]);set(gca,'yTickLabel',[tx,ty,1],fss,afs)
title('rayleigh score per cell'); 
xlabel('pre',fss,afs); ylim([0 1]); xlim([0 1]);
ylabel('dur',fss,afs); axis square; hold off;

set(gC,'Heights', [h -1]);

%% D pairs time by hd
gD = uix.GridFlex(p,gCent,spss,1, bcss,'w',dmss,ofss);
axes(p,uicontainer(p,gD,bcss,'w'),vss,ofss);
%Text
text(0,0.9,'D',fwss,bss,fss,fs);
axes(p,uicontainer(p,gD,bcss,'w'),vss,ofss);
tl1=cl1&rmi; tl2=cl2&rmi;
plot(ctsbma(tl1,1), ctsbma(tl1,2),'b.',mss,mfs,'linewidth',lfs); hold on; %MISSING SOME??
plot(ctsbma(tl2,1), ctsbma(tl2,2),'r.',mss,mfs); 
f1 = fit(ctsbma(tl1,1),ctsbma(tl1,2),'poly1');
mm = slimd(gca); pscat = plot(mm,f1(mm),'b:'); pscat.LineWidth = lfs-1;
f2 = fit(ctsbma(tl2,1),ctsbma(tl2,2),'poly1');
mm = slimd(gca); pscat = plot(mm,f2(mm),'r:'); pscat.LineWidth = lfs-1;
text(0.1,0.9,...
sprintf('r=%.2f a=%.2f',rnd(ccof(ctsbma(tl1,1),ctsbma(tl1,2)),2),rnd(f1.p1,2)),coss,'b',unss,nzss);
text(0.1,0.8,...
sprintf('r=%.2f a=%.2f',rnd(ccof(ctsbma(tl2,1),ctsbma(tl2,2)),2),rnd(f2.p1,2)),coss,'r',unss,nzss);
title('temporal correlations',fss,tfs); legend off;
xlabel('pre',fss,afs); ylabel('dur',fss,afs);xlim(mm); ylim(mm); axis square;
set(gD,'Heights', [h -1]);

%% E pairs spatial by hd
gE = uix.GridFlex(p,gCent,spss,1, bcss,'w',dmss,ofss);
axes(p,uicontainer(p,gE,bcss,'w'),vss,ofss);
text(0,0.9,'E',fwss,bss,fss,fs);
axes(p,uicontainer(p,gE,bcss,'w'),vss,ofss);
%plot(ctsbma(rmi,4),   ctsbma(rmi,5),'kx',mss,mfs,'linewidth',lfs-.5);hold on;
tl1=cl1&rmi; tl2=cl2&rmi;
plot(ctsbma(tl1,4), ctsbma(tl1,5),'b.',mss,mfs,'linewidth',lfs);hold on;
plot(ctsbma(tl2,4), ctsbma(tl2,5),'r.',mss,mfs);
f1 = fit(ctsbma(tl1,4),ctsbma(tl1,5),'poly1');
mm = slimd(gca); pscat = plot(mm,f1(mm),'b:'); pscat.LineWidth = lfs-1;
f2 = fit(ctsbma(tl2,4),ctsbma(tl2,5),'poly1');
mm = slimd(gca); pscat = plot(mm,f2(mm),'r:'); pscat.LineWidth = lfs-1;
text(0.1,0.9,...
sprintf('r=%.2f a=%.2f',round(ccof(ctsbma(tl1,4),ctsbma(tl1,5)),2),round(f1.p1,2)),coss,'b',unss,nzss);
text(0.1,0.8,...
sprintf('r=%.2f a=%.2f',round(ccof(ctsbma(tl2,4),ctsbma(tl2,5)),2),round(f2.p1,2)),coss,'r',unss,nzss);
title('spatial correlations',fss,tfs); legend off;
xlabel('pre',fss,afs);ylabel('dur',fss,afs);xlim(mm); ylim(mm); axis square;
set(gE,'Heights', [h -1]);

%% gCent end
set(gCent,'Heights', [-1 -1 -1 -1]);

%% Right
gRight = uix.GridFlex(p, gTop,spss,5, bcss,'w',dmss,ofss);

%% F hd cell 
gF = uix.GridFlex(p,gRight,spss,1, bcss,'w',dmss,ofss);
axes(p,uicontainer(p,gF,bcss,'w'),vss,ofss);
text(tx,ty2,'F',fwss,bss,fss,fs);
rst = 0.3; yl=12;
axes(p,uicontainer(p,gF,bcss,'w'),vss,ofss); 
% title(sprintf('r-angle hd cells r-score > %0.1f',rst));
x = 2; t = []; a=fp(x, hdclls, rst,yl,'hd cluster dur','r'); len(t);
set(gF,'Heights', [h -1]);    

        
%% G hd cell by session
gG = uix.GridFlex(p,gRight,spss,1, bcss,'w',dmss,ofss);
axes(p,uicontainer(p,gG,bcss,'w'),vss,ofss);
text(tx,ty2,'G',fwss,bss,fss,fs);
rst = 0.3; t1 = 'r-angle '; t2=' all';yl=41;
%title(sprintf('r-angle for all cells with r-score > %0.1f',rst));
axes(p,uicontainer(p,gG,bcss,'w'),vss,ofss);
x = 1; t = [];fp(x,cellsn,rst,yl,[t1 tss{x} t2]); len(t); 
uicontainer(p,gG,bcss,'w'); %empty;
axes(p,uicontainer(p,gG,bcss,'w'),vss,ofss);  
x = 2; t = [];fp(x,cellsn,rst,yl,[t1 tss{x} t2]); len(t);
uicontainer(p,gG,bcss,'w'); %empty;
axes(p,uicontainer(p,gG,bcss,'w'),vss,ofss);  
x = 3; t = [];fp(x,cellsn,rst,yl,[t1 tss{x} t2]); len(t);        
function h=fp(x,cll,rst,yl,tit,c)
    cla;
    t = [cll.(s{x})]; tl=arrayfun(@(z) len(z.st),t); tl=tl>=100;   
    ta = rad2deg([t.(va)]); ts = [t.(vs)];  tis = ts>rst; t = ta(tis & tl);
    h=histogram(t,-180:30:180); xlim([-180 180]);  axis square; ylim([0 yl]);
    h.FaceColor='g'; if nargin==6; h.FaceColor=c; h.FaceAlpha=1; end 
    title(tit,fss,tfs); xlabel('r-angle'); ylabel(' ');
    l=[-180 180]; d=char(176); sf='%d%c';
    set(gca,'xtick',l,'xticklabel',{sprintf(sf,l(1),d); sprintf(sf,l(2),d)});
    text(0.1,0.95,sprintf('n=%d a=%.2f',len(t),meanang(t)),unss,nzss);
end        
set(gG,'Heights', [h -1 h -1 h -1]);    


%% gRight
set(gRight,'Heights', [-1 -3]);




%% epilogue
set(gTop,'width', [-3 -1 -1]);

end
%
