% uncomment to load cells
%load('Z:\\data\\noam\\Data\\muscimol\\noam\\cells_45min_d_patchtraj_rayleigh'); %lab
%load('C:\\Noam\\Data\\muscimol\\noam\\cells_45min_d_patchtraj_rayleigh'); %personal

%%preamble 
dbstop if error  
%loads();
f = figure(991);
set(f,'Color','w', 'Position', [600 0 1200 800]);
gtop = uix.GridFlex('Parent',f,'Spacing',5, 'BackgroundColor','w');
gl = uix.GridFlex('Parent',gtop,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
axes('Parent',uicontainer('Parent',gl,'BackgroundColor','w'),'visible','off');
text(0,0.5,'A','fontweight','bold','fontsize',fs);
% uix.Empty('Parent',gl);
% txt = uicontrol('Style','text', 'String','A','Parent',gl,...
%     'BackgroundColor','w','FontSize',11, 'fontweight','bold','position',[0 0 0 0])
gA = uix.GridFlex('Parent',gl,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
%axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w'),'visible','off');
%text(0,0.5,'A','fontweight','bold','fontsize',fs);

set(gl,'Widths', [-1], 'Heights', [ 25, -1]);
gr = uix.GridFlex('Parent',gtop,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');

%1A - trajectory / auto for group
ii = [100, 103, 104, 105, 111];
axg = []; ug = []; s1 = []; s2 = s1;
s1.x='position(cm)'; s1.y = s1.x;
s2.x=''; s2.y = s2.x;
for i = 1:len(ii)
    c = cells{ii(i)};
    s1.t = ' ';
    axg(end+1) = axes('Parent',uicontainer('Parent',gA)); if(i==3); s1.t = 'pre'; end 
    axg(end)=plotTR(axg(end),c.before,s1);
    s1.x=''; s1.y = s1.x;
    axg(end+1) = axes('Parent',uicontainer('Parent',gA));  if(i==3); s1.t = 'dur'; end 
    axg(end)=plotTR(axg(end),c.midall,s1);
    axg(end+1) = axes('Parent',uicontainer('Parent',gA));  if(i==3); s1.t = 'post'; end 
    axg(end)=plotTR(axg(end),c.after,s1);
    axg(end+1) = axes('Parent',uicontainer('Parent',gA));  s2.t = sprintf('%s%.1f','grid=',c.before.gridscore);
    axg(end)=plotAC(axg(end),c.before,s2);
    axg(end+1) = axes('Parent',uicontainer('Parent',gA));  s2.t = sprintf('%s%.1f','grid=',c.midall.gridscore);
    axg(end)=plotAC(axg(end),c.midall,s2);
    axg(end+1) = axes('Parent',uicontainer('Parent',gA));  s2.t = sprintf('%s%.1f','grid=',c.after.gridscore);
    axg(end)=plotAC(axg(end),c.after,s2);
end
set(gA,'Widths', [-1 -1 -1 -1 -1], 'Heights', [-1 -1 -1 -1 -1 -1])

%1B - mean firing pre dur
axes('Parent',uicontainer('Parent',gr,'BackgroundColor','w'),'visible','off');
text(0,0.5,'B','fontweight','bold','fontsize',fs);
uB = uicontainer('Parent',gr); axB = axes(uB); mrbda = [];
for i = 1:len(cels)
    c = cells{cels(i)};
    mrbda(i,1) = len(c.before.st)/(c.before.st(end)-c.before.st(1));
    mrbda(i,2) = len(c.midall.st)/(c.midall.st(end)-c.midall.st(1));
    mrbda(i,3) = len(c.after.st)/(c.after.st(end)-c.after.st(1));
end
plot(mrbda(:,1),mrbda(:,2),'ko'); hold on;
% plot(mrbda(:,3),mrbda(:,2),'bo');
hold off
axB.XLim = slim(axB); axB.YLim = axB.XLim; axis(axB,'square');
xlabel('pre (Hz)'); ylabel('dur (Hz)');
title('mean firing rate pre during');
%legend(axB, {'pre','post'});


%1D threshold all
axes('Parent',uicontainer('Parent',gr,'BackgroundColor','w'),'visible','off');
text(0,0.5,'D','fontweight','bold','fontsize',fs);
uD = uicontainer('Parent',gr); axD = axes(uD);
sigbma = [];
for i = 1:len(cels)
    c = cells{cels(i)};
    sigbma(i,:) = [c.before.gridscore c.midall.gridscore c.after.gridscore];     
end
agbm = [];
for i = 1:len(cells)
    c = cells{i};
    agbm(i,:) = [c.before.gridscore c.midall.gridscore];     
end
%convert to 0
agbm(agbm == -2) = 0;
sigbma(sigbma == -2) = 0;
plot(agbm(:,1),agbm(:,2),'b.'); hold on
plot(sigbma(:,1),sigbma(:,2),'ro');
%plot(gbm((gbm(:,2) == 0),1),gbm((gbm(:,2) == 0),2),'rx');
axD.XLim = slimd(axD); axD.YLim = axD.XLim;
plot([0.3 0.3],[-1 0.25],'g');plot([0.3 2],[0.25 0.25],'g');
hold off; axis(axD,'square');
xlabel('gridscore pre'); ylabel('gridscore during');
title('gridscore all vs cohort'); 
%legend(axC, {'all','significant','thresholds'});

%1C - mean firing pre post
axes('Parent',uicontainer('Parent',gr,'BackgroundColor','w'),'visible','off');
text(0,0.5,'C','fontweight','bold','fontsize',fs);
uC = uicontainer('Parent',gr); axC = axes(uC); 
plot(mrbda(:,1),mrbda(:,3),'ko'); hold on;
hold off
axC.XLim = slim(axC); axC.YLim = axC.XLim; axis(axC,'square');
xlabel('pre (Hz)'); ylabel('post (Hz)');
title('mean firing rate pre post');

%1E Significance cels only
axes('Parent',uicontainer('Parent',gr,'BackgroundColor','w'),'visible','off');
text(0,0.5,'E','fontweight','bold','fontsize',fs);
uE = uicontainer('Parent',gr,'BackgroundColor','w');axE = axes(uE); 
%figure(9914)aspace
plot(sigbma(:,1),sigbma(:,3),'ko');
hold off
axE.XLim = slim(axE); axE.YLim = axE.XLim; axis(axE,'square');
xlabel('gridscore pre'); ylabel('gridscore post');
title('gridscore pre vs post');
set( gr, 'Heights', [25 -1 25 -1] );
%set( gr, 'Widths', [-1], 'Heights', [-1 -1 -1]);
a = findobj(f,'type','UIContainer');
for i = 1:len(a)
    a(i).BackgroundColor = 'w';
end

%print(f, 'f1.png', '-dpng', '-r300' );   %save file as PNG w/ 300dpi

