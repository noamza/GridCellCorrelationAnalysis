% uncomment to load cells
%load(sprintf('C:\\Noam\\Data\\muscimol\\noam\\cells_%dmin_d_patchtraj_rayleigh',45));
%load('Z:\\data\\noam\\Data\\muscimol\\noam\\cells_45min_d_patchtraj_rayleigh'); %lab
%%{
load('.\data\pairs'); cels = unique(pairs(:))';
load('.\data\cptsbma'); rmi = []; %rmi=[43];

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
ccl2 = rs(:,2)>0.4; ccl1 = ~ccl2; chd = cels(ccl2);
cl2 = rsp(:,3)>0.4 & rsp(:,4)>0.4; cl1 = ~cl2;

%}

%preamble
dbstop if error    
f = figure(99433); tic
set(f,'Color','w', 'Position', [200 0 800 900]); pf = {};
gtop = uix.GridFlex('Parent', f,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
gleft = uix.GridFlex('Parent', gtop,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
ci = [112, 114, 117, 118, 119]; %ci = [103, 104, 105, 100,111];  ci = [282,283,287]; ci = [112, 114, 115, 116, 117,118,119]; 115
fs = 12; afs = 10; mfs = 10; tfs = 12; lfs = 3; h = 25;

%% A
%axes('Parent',uicontainer('Parent',gleft,'BackgroundColor','w'),'visible','off');
axes('Parent',uix.Panel('Parent',gleft,'BackgroundColor','w','BorderType','none'),'visible','off');
text(0,0.5,'A','fontweight','bold','fontsize',fs); axis off;
gA = uix.GridFlex('Parent',gleft,'Spacing',1, 'BackgroundColor','w','DividerMarkings','on');
axA = {}; pf = []; s1 = [];
for i = 1:len(ci)
    c = cells{ci(i)};
    %ggA = uix.VBox('Parent',gA,'Spacing',1, 'BackgroundColor','w');%,'DividerMarkings','off');
    %TRAJ
    axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w'),'visible','off');
    text(0.5,0.5,sprintf('C%d',i),'fontweight','bold','fontsize',mfs,'HorizontalAlignment','center');
    s1.x=''; s1.y = s1.x; s1.t ='';
    axA{end+1} = axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w'));
    axA{end} = plotTR(axA{end},c.before,s1); s1.t = '';
    axA{end+1} = axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w'));
    axA{end} =plotTR(axA{end},c.midall,s1);
    %AC
    %uix.Empty('Parent',ggA); 
    s1.t = sprintf('grid=%.1f',round(c.before.gridscore,2));
    axA{end+1} = axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w'));
    axA{end} = plotAC(axA{end},c.before,s1);
    s1.t = sprintf('grid=%.1f',round(c.midall.gridscore,2));
    axA{end+1} = axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w'));
    axA{end} = plotAC(axA{end},c.midall,s1);
    %RAY
    %uix.Empty('Parent',ggA); 
    s1.t = sprintf('rayleigh=%.1f',round(c.before.rayleigh_score,2));
    axA{end+1} = axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w'));
    axA{end} = plotHD(axA{end},c.before,s1);
    s1.t = sprintf('rayleigh=%.1f',round(c.midall.rayleigh_score,2));
    axA{end+1} = axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w'));
    axA{end} = plotHD(axA{end},c.midall,s1);
end
%axA{end+1} = axes('Parent',uicontainer('Parent',ggA,'BackgroundColor','w'),'visible','off');
%uix.Empty('Parent',ggA);
set(gA,'Heights', [15 -1 -1 -1 -1 -1 -1]);%,'Widths', [-1 -1 -1]);


%% B
axes('Parent',uix.Panel('Parent',gleft,'BackgroundColor','w','BorderType','none'),'visible','off');
text(0,0.0,'B','fontweight','bold','fontsize',fs,'HorizontalAlignment','center');  axis off;
gB = uix.GridFlex('Parent',gleft,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
axB = {};
s1.x=''; s1.y = ''; s1.t = ''; pf.off = 0; pf.lag = 1000; pf.movmean = 100;
%%{
for i = 1:len(ci)-1
    for ii = i+1:len(ci)
        c1 = cells{ci(i)};c2 = cells{ci(ii)};
        if i+ii>3
            %uix.Empty('Parent',gB);
        end
        %s1.t = sprintf('C%dxC%d',i,ii);
        axes('Parent',uicontainer('Parent',gB,'BackgroundColor','w'),'visible','off');
        text(0.3,0.5,sprintf('C%dxC%d',i,ii),'fontweight','bold','fontsize',mfs,'HorizontalAlignment','center');
        axB{end+1} = axes('Parent',uicontainer('Parent',gB,'BackgroundColor','w'));
        axB{end}  = plotTimeCorr(axB{end},c1.before,c2.before,pf,s1);
        axB{end+1} = axes('Parent',uicontainer('Parent',gB,'BackgroundColor','w'));
        s1.t = '';
        axB{end}  = plotTimeCorr(axB{end},c1.midall,c2.midall,pf,s1); 
        axB{end+1} = axes('Parent',uicontainer('Parent',gB,'BackgroundColor','w')); 
        axB{end}  = plotSpaceCorr(axB{end},c1.before,c2.before,pf,s1);
        axB{end+1} = axes('Parent',uicontainer('Parent',gB,'BackgroundColor','w')); 
        axB{end}  = plotSpaceCorr(axB{end},c1.midall,c2.midall,pf,s1);
        
    end
end
set(gB,'Heights', [15 -1 -1 -1 -1]);%,'Widths', [-1 -1 -1]);
%}
%gleft
set(gleft,'Heights', [25 -2 25 -1]);

%gright
gright = uix.GridFlex('Parent', gtop,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
%% C hd single
gA = uix.GridFlex('Parent',gright,'Spacing',1, 'BackgroundColor','w','DividerMarkings','off');
axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w'),'visible','off');
%Text
text(0,0.9,'C','fontweight','bold','fontsize',fs);
axA = axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w'),'visible','off');
plot(rs(:,1),rs(:,2),'b.','markersize',mfs);
%hold on; plot(rs(rs(:,2)>0.5,1),rs(rs(:,2)>0.5,2),'r.','markersize',mfs);
hold on; plot(rs(ccl2,1),rs(ccl2,2),'r.','markersize',mfs);

%title('rayleigh score by cell');
legend({'cluster 1';'cluster 2'},'fontsize',afs,'Location','east');
set(gca,'XTick',[0,0.5,1]);set(gca,'XTickLabel',[0,0.5,1],'fontsize',afs)
set(gca,'YTick',[0,0.5,1]);set(gca,'yTickLabel',[0,0.5,1],'fontsize',afs)
title('rayleigh score per cell');
xlabel('pre','fontsize',afs); axA.YLim=[0 1];axA.XLim=[0 1];
ylabel('dur','fontsize',afs);axis(axA,'square');
set(gA,'Heights', [h -1]);
%% D time by hd
gB = uix.GridFlex('Parent',gright,'Spacing',1, 'BackgroundColor','w','DividerMarkings','off');
axes('Parent',uicontainer('Parent',gB,'BackgroundColor','w'),'visible','off');
%Text
text(0,0.9,'D','fontweight','bold','fontsize',fs);
axB = axes('Parent',uicontainer('Parent',gB,'BackgroundColor','w'),'visible','off');
%plot(ctsbma(rmi,1),   ctsbma(rmi,2),'kx','markersize',mfs,'linewidth',lfs-.5);hold on;  %REMOVE OUTLIER
plot(ctsbma(:,1),   ctsbma(:,2),'b.','markersize',mfs,'linewidth',lfs); hold on;
plot(ctsbma(cl2,1), ctsbma(cl2,2),'r.','markersize',mfs);
tc1=cl1;  tc2=cl2; %tc1(rmi)=false; tc2(rmi)=false; %REMOVE OUTLIER
f1 = fit(ctsbma(tc1,1),ctsbma(tc1,2),'poly1');
mm = slimd(gca); pscat = plot(mm,f1(mm),'b:'); pscat.LineWidth = lfs-1;
f2 = fit(ctsbma(tc2,1),ctsbma(tc2,2),'poly1');
mm = slimd(gca); pscat = plot(mm,f2(mm),'r:'); pscat.LineWidth = lfs-1;
text(0.1,0.9,sprintf('r=%.2f a=%.2f',round(ccof(ctsbma(tc1,1),ctsbma(tc1,2)),2),round(f1.p1,2)),...
    'color','b','Units','normalized');
text(0.1,0.8,sprintf('r=%.2f a=%.2f',round(ccof(ctsbma(tc2,1),ctsbma(tc2,2)),2),round(f2.p1,2)),...
    'color','r','Units','normalized');
%text(0.1,0.7,sprintf('r=%.2f',round(ccof(ctsbma(pptsbma(:,2),1),ctsbma(pptsbma(:,2),2)),2)),'color','m','Units','normalized');
xlabel('pre','fontsize',afs);
ylabel('dur','fontsize',afs)
title('temporal correlations','fontsize',tfs); legend off;
axB.XLim = mm; axB.YLim = axB.XLim; axis(axB,'square');
set(gB,'Heights', [h -1]);
%% E spatial by hd
gC = uix.GridFlex('Parent',gright,'Spacing',1, 'BackgroundColor','w','DividerMarkings','off');
axes('Parent',uicontainer('Parent',gC,'BackgroundColor','w'),'visible','off');
%Text
text(0,0.9,'E','fontweight','bold','fontsize',fs);
axC = axes('Parent',uicontainer('Parent',gC,'BackgroundColor','w'),'visible','off');
%plot(ctsbma(rmi,4),   ctsbma(rmi,5),'kx','markersize',mfs,'linewidth',lfs-.5);hold on;
plot(ctsbma(:,4),   ctsbma(:,5),'b.','markersize',mfs,'linewidth',lfs);hold on;
plot(ctsbma(cl2,4), ctsbma(cl2,5),'r.','markersize',mfs);
tc1=cl1; tc1(rmi)=false; tc2=cl2; tc2(rmi)=false;
f1 = fit(ctsbma(cl1,4),ctsbma(cl1,5),'poly1');
mm = slimd(gca); pscat = plot(mm,f1(mm),'b:'); pscat.LineWidth = lfs-1;
f2 = fit(ctsbma(cl2,4),ctsbma(cl2,5),'poly1');
mm = slimd(gca); pscat = plot(mm,f2(mm),'r:'); pscat.LineWidth = lfs-1;
text(0.1,0.9,sprintf('r=%.2f a=%.2f',round(ccof(ctsbma(cl1,4),ctsbma(cl1,5)),2),round(f1.p1,2)),...
    'color','b','Units','normalized');
text(0.1,0.8,sprintf('r=%.2f a=%.2f',round(ccof(ctsbma(cl2,4),ctsbma(cl2,5)),2),round(f2.p1,2)),...
    'color','r','Units','normalized');
%text(0.1,0.7,sprintf('r=%.2f',round(ccof(ctsbma(pptsbma(:,5),4),ctsbma(pptsbma(:,5),5)),2)),'color','m','Units','normalized');
xlabel('pre','fontsize',afs);
ylabel('dur','fontsize',afs)
title('spatial correlations','fontsize',tfs); legend off;
axC.XLim = mm; axC.YLim = axC.XLim; axis(axC,'square');
set(gC,'Heights', [h -1]);


%% gright
set(gright,'Heights', [-1 -1 -1]);
%% epilogue
set(gtop,'width', [-3 -1]);

%
