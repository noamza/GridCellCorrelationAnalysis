% uncomment to load cells
%load(sprintf('C:\\Noam\\Data\\muscimol\\noam\\cells_%dmin_d_patchtraj_rayleigh',45));
%load('Z:\\data\\noam\\Data\\muscimol\\noam\\cells_45min_d_patchtraj_rayleigh'); %lab

%% preamble
dbstop if error
%%loads();
%%pairs(43,:)=[];
%{
rsp=[];rs = [];
for i=1:length(pairs)
    rsp(i,:) = [cells{pairs(i,1)}.before.rayleigh_score,...
        cells{pairs(i,2)}.before.rayleigh_score,...
        cells{pairs(i,1)}.midall.rayleigh_score,...
        cells{pairs(i,2)}.midall.rayleigh_score];
end
for i=1:length(cels)
    rs(i,:)  = [cells{cels(i)}.before.rayleigh_score,...
        cells{cels(i)}.midall.rayleigh_score];
end
cl2 = rsp(:,3)>0.4 & rsp(:,4)>0.4; cl1 = ~cl2;

% pairs(43,:)=[];ctsbma(43,:)=[]; ptsbma(43,:)=[];pptsbma(43,:)=[];cl1(43)=[];cl2(43)=[];zpptsbma(43,:)=[];
%}
f = figure(994); tic
set(f,'Color','w', 'Position', [200 0 1200 900]);
fs = 12; afs = 10; mfs = 10; tfs = 12; lfs = 3; p = {};
gtop = uix.GridFlex('Parent', f,'Spacing',5, 'BackgroundColor','w','DividerMarkings','on');

% %% A
% gA = uix.GridFlex('Parent',gtop,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
% axA = {};
% axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w'),'visible','off');
% %Text
% text(0,0.5,'A','fontweight','bold','fontsize',fs);
% axA{end+1} = axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w'),'visible','off');
% text(0.5,0.5,'pre','fontsize',fs-1,'fontweight','bold');
% axA{end+1} = axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w'),'visible','off');
% text(0.5,0.5,'dur','fontsize',fs-1,'fontweight','bold');
% %TRAJ 1
% uix.Empty('Parent',gA); i = 42;
% c1 = cells{pairs(i,1)}; c2 = cells{pairs(i,2)};
% s1.x=''; s1.y = s1.x; s1.t = sprintf('C%d',c1.ind);
% axA{end+1} = axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w'));
% axA{end} = plotTR(axA{end},c1.before,s1); s1.t = '';
% axA{end+1} = axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w'));
% axA{end} =plotTR(axA{end},c1.midall,s1);
% %AC 1
% uix.Empty('Parent',gA); 
% s1.t = sprintf('grid=%.1f',round(c1.before.gridscore,2));
% axA{end+1} = axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w'));
% axA{end} = plotAC(axA{end},c1.before,s1);
% s1.t = sprintf('grid=%.1f',round(c1.midall.gridscore,2));
% axA{end+1} = axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w'));
% axA{end} = plotAC(axA{end},c1.midall,s1);
% %ray 1
% uix.Empty('Parent',gA); 
% s1.t = sprintf('rayleigh=%.1f',round(c1.before.rayleigh_score,2));
% axA{end+1} = axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w'));
% axA{end} = plotHD(axA{end},c1.before,s1);
% s1.t = sprintf('rayleigh=%.1f',round(c1.midall.rayleigh_score,2));
% axA{end+1} = axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w'));
% axA{end} = plotHD(axA{end},c1.midall,s1);
% %TRAJ 2
% uix.Empty('Parent',gA);s1.x=''; s1.y = s1.x; s1.t = sprintf('C%d',c2.ind);
% axA{end+1} = axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w'));
% axA{end} = plotTR(axA{end},c2.before,s1); s1.t ='';
% axA{end+1} = axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w'));
% axA{end} =plotTR(axA{end},c2.midall,s1);
% %AC 2
% uix.Empty('Parent',gA); 
% s1.t = sprintf('grid=%.1f',round(c2.before.gridscore,2));
% axA{end+1} = axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w'));
% axA{end} = plotAC(axA{end},c2.before,s1);
% s1.t = sprintf('grid=%.1f',round(c2.midall.gridscore,2));
% axA{end+1} = axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w'));
% axA{end} = plotAC(axA{end},c2.midall,s1);
% %ray 2
% uix.Empty('Parent',gA);
% axA{end+1} = axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w'));
% s1.t = sprintf('rayleigh=%.1f',round(c2.before.rayleigh_score,2));
% axA{end} = plotHD(axA{end},c2.before,s1);
% axA{end+1} = axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w'));
% s1.t = sprintf('rayleigh=%.1f',round(c2.midall.rayleigh_score,2));
% axA{end} = plotHD(axA{end},c2.midall,s1);
% %TCC
% p.lag = 2000; p.movmean = 100; p.ylim = [-.06 0.06];
% uix.Empty('Parent',gA); s1.t = 'xcorr';
% axA{end+1} = axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w'));
% axA{end}  = plotTimeCorr(axA{end},c1.before,c2.before,p,s1);
% axA{end+1} = axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w'));
% axA{end}  = plotTimeCorr(axA{end},c1.midall,c2.midall,p,s1);
% %uix.Empty('Parent',gA);
% set(gA,'Heights', [25 -1 -1]);%,'Widths', [-1 -1 -1]);

%%
gmid = uix.GridFlex('Parent', gtop,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off'); 
w=25;

%% B hd pairs pre
gB = uix.GridFlex('Parent',gmid,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
axes('Parent',uicontainer('Parent',gB,'BackgroundColor','w'),'visible','off');
%Text
text(0,0.9,'B','fontweight','bold','fontsize',fs);
axB = axes('Parent',uicontainer('Parent',gB,'BackgroundColor','w'),'visible','off');
plot(rsp(:,1), rsp(:,2),'k.','markersize',mfs); %plot(ax(end),d.x2, d.y2, 'ro','linewidth',1.2);
title('head direction pre','FontSize',tfs);
xlabel('cell 1','FontSize',afs);ylabel('cell 2','FontSize',afs);
xlim([0 1]); ylim([0 1]); axis 'square';
set(gB,'Widths', [w -1]);

%% E time by hd
gE = uix.GridFlex('Parent',gmid,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
axes('Parent',uicontainer('Parent',gE,'BackgroundColor','w'),'visible','off');
%Text
text(0,0.9,'E','fontweight','bold','fontsize',fs);
axE = axes('Parent',uicontainer('Parent',gE,'BackgroundColor','w'),'visible','off');
plot(ctsbma(:,1),   ctsbma(:,2),'b.','markersize',mfs,'linewidth',lfs); hold on;
plot(ctsbma(cl2,1), ctsbma(cl2,2),'r.','markersize',mfs);
f1 = fit(ctsbma(cl1,1),ctsbma(cl1,2),'poly1');
mm = slimd(gca); pscat = plot(mm,f1(mm),'b:'); pscat.LineWidth = lfs-1;
f2 = fit(ctsbma(cl2,1),ctsbma(cl2,2),'poly1');
mm = slimd(gca); pscat = plot(mm,f2(mm),'r:'); pscat.LineWidth = lfs-1;
%plot(ctsbma(pptsbma(:,2),1), ctsbma(pptsbma(:,2),2),'mo','markersize',mfs,'linewidth',lfs-2.5);
%plot(ctsbma(pptsbma(:,1),1), ctsbma(pptsbma(:,1),2),'go','markersize',mfs,'linewidth',lfs-2.5);
%f3 = fit(ctsbma(pptsbma(:,2),1), ctsbma(pptsbma(:,2),2),'poly1');
%textbest(axC,sprintf('%.2fx + %.2f',round(f1.p1,2),round(f1.p2,2)),{},{}); 'a=%.2fx+%.2f'
% text(0.1,0.9,sprintf('%.2fx + %.2f',round(f1.p1,2),round(f1.p2,2)),'color','b','Units','normalized');
% text(0.1,0.8,sprintf('%.2fx + %.2f',round(f2.p1,2),round(f2.p2,2)),'color','r','Units','normalized');
% text(0.1,0.7,sprintf('%.2fx + %.2f',round(f3.p1,2),round(f3.p2,2)),'color','m','Units','normalized');
text(0.1,0.9,sprintf('r=%.2f',round(ccof(ctsbma(cl1,1),ctsbma(cl1,2)),2)),'color','b','Units','normalized');
text(0.1,0.8,sprintf('r=%.2f',round(ccof(ctsbma(cl2,1),ctsbma(cl2,2)),2)),'color','r','Units','normalized');
%text(0.1,0.7,sprintf('r=%.2f',round(ccof(ctsbma(pptsbma(:,2),1),ctsbma(pptsbma(:,2),2)),2)),'color','m','Units','normalized');
xlabel('pre','fontsize',afs);
ylabel('dur','fontsize',afs)
title('time correlations','fontsize',tfs); legend off;
axE.XLim = mm; axE.YLim = axE.XLim; axis(axE,'square');
set(gE,'Widths', [w -1]);


%% C hd pairs dur
gC = uix.GridFlex('Parent',gmid,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
axes('Parent',uicontainer('Parent',gC,'BackgroundColor','w'),'visible','off');
%Text
text(0,0.9,'C','fontweight','bold','fontsize',fs);
axC = axes('Parent',uicontainer('Parent',gC,'BackgroundColor','w'),'visible','off');
plot(rsp(:,3), rsp(:,4),'k.','markersize',mfs);
title('head direction dur','FontSize',tfs);
xlabel('cell 1','FontSize',afs);ylabel('cell 2','FontSize',afs)
hold on;
plot([0.4 0.4 1],[1 0.4 0.4],'r--','linewidth',lfs); %plot(ax(end),[0.51 0.51 1],[1 0.51 0.51],'r--');
text(gca,0.1,0.95,'cluster 1','Color','b','FontSize',afs);
text(gca,0.6,0.95,'cluster 2','Color','r','FontSize',afs);
xlim([0 1]); ylim([0 1]); axis 'square';
set(gC,'Widths', [w -1]);


%% F spatial by hd
gF = uix.GridFlex('Parent',gmid,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
axes('Parent',uicontainer('Parent',gF,'BackgroundColor','w'),'visible','off');
%Text
text(0,0.9,'F','fontweight','bold','fontsize',fs);
axF = axes('Parent',uicontainer('Parent',gF,'BackgroundColor','w'),'visible','off');
plot(ctsbma(:,4),   ctsbma(:,5),'b.','markersize',mfs,'linewidth',lfs); hold on;
plot(ctsbma(cl2,4), ctsbma(cl2,5),'r.','markersize',mfs);
f1 = fit(ctsbma(cl1,4),ctsbma(cl1,5),'poly1');
mm = slimd(gca); pscat = plot(mm,f1(mm),'b:'); pscat.LineWidth = lfs-1;
f2 = fit(ctsbma(cl2,4),ctsbma(cl2,5),'poly1');
mm = slimd(gca); pscat = plot(mm,f2(mm),'r:'); pscat.LineWidth = lfs-1;
%plot(ctsbma(pptsbma(:,5),4), ctsbma(pptsbma(:,5),5),'mo','markersize',mfs,'linewidth',lfs-2.5);
%f3 = fit(ctsbma(pptsbma(:,5),4), ctsbma(pptsbma(:,5),5),'poly1');
%textbest(axC,sprintf('%.2fx + %.2f',round(f1.p1,2),round(f1.p2,2)),{},{}); 'a=%.2fx+%.2f'
% text(0.1,0.9,sprintf('%.2fx + %.2f',round(f1.p1,2),round(f1.p2,2)),'color','b','Units','normalized');
% text(0.1,0.8,sprintf('%.2fx + %.2f',round(f2.p1,2),round(f2.p2,2)),'color','r','Units','normalized');
% text(0.1,0.7,sprintf('%.2fx + %.2f',round(f3.p1,2),round(f3.p2,2)),'color','m','Units','normalized');
text(0.1,0.9,sprintf('r=%.2f',round(ccof(ctsbma(cl1,4),ctsbma(cl1,5)),2)),'color','b','Units','normalized');
text(0.1,0.8,sprintf('r=%.2f',round(ccof(ctsbma(cl2,4),ctsbma(cl2,5)),2)),'color','r','Units','normalized');
%text(0.1,0.7,sprintf('r=%.2f',round(ccof(ctsbma(pptsbma(:,5),4),ctsbma(pptsbma(:,5),5)),2)),'color','m','Units','normalized');
xlabel('pre','fontsize',afs);
ylabel('dur','fontsize',afs)
title('spatial correlations','fontsize',tfs); legend off;
axF.XLim = mm; axF.YLim = axF.XLim; axis(axF,'square');
set(gF,'Widths', [w -1]);

%% D hd single
gD = uix.GridFlex('Parent',gmid,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
axes('Parent',uicontainer('Parent',gD,'BackgroundColor','w'),'visible','off');
%Text
text(0,0.9,'D','fontweight','bold','fontsize',fs);
axD = axes('Parent',uicontainer('Parent',gD,'BackgroundColor','w'),'visible','off');
plot(rs(:,1),rs(:,2),'b.','markersize',mfs);
hold on; plot(rs(rs(:,2)>0.5,1),rs(rs(:,2)>0.5,2),'r.','markersize',mfs);
%title('rayleigh score by cell');
legend({'cluster 1';'cluster 2'},'fontsize',afs,'Location','best');
set(gca,'XTick',[0,0.5,1]);set(gca,'XTickLabel',[0,0.5,1],'fontsize',afs)
set(gca,'YTick',[0,0.5,1]);set(gca,'yTickLabel',[0,0.5,1],'fontsize',afs)
title('rayleigh score per cell');
xlabel('pre','fontsize',afs);
ylabel('dur','fontsize',afs);axis(axD,'square');
set(gD,'Widths', [w -1]);

set(gmid,'Heights', [-1, -1],'Widths', [-1 -1 -1]);



%% epilogue
set(gtop,'Heights', [-1]);%,'Widths', [ -1 -1]);

%}