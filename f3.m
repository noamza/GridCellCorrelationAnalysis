%% preamble
dbstop if error  
%load('.\\data\\shuffling1000nan','ctsbma','ptsbma','pptsbma');
nshuffle=1000;
fig = figure(993); tic
set(fig,'Color','w', 'Position', [200 0 1000 600]);
fss='fontsize'; mss='markersize';fw={'fontweight','bold'}; pnt='Paren';
fs ={fss 14}; afs ={fss 10}; mfs ={mss 10}; tfs = {fss,12}; 


gTop = uix.GridFlex('Parent', fig,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
gfp={'Spacing',5, 'BackgroundColor','w','DividerMarkings','off'};

tx=0;ty=0.5;tp={fw{:},fs{:}};
uip={'BackgroundColor','w'};
aup={'visible','off'};
%{
'sig:sig   sig:non   non:sig   non:non'
[sum(pptsbma(:,a)&pptsbma(:,b)),...
sum(pptsbma(:,a)&~pptsbma(:,b)),...
sum(~pptsbma(:,a)&pptsbma(:,b)),...
sum(~pptsbma(:,a)&~pptsbma(:,b)) ]
if(a==1) disp('pre'); else disp('dur'); end
%}

%% A1 : Before During Temporal
gA1 = uix.GridFlex('Parent',gTop,gfp{:});
axes('Parent',uicontainer('Parent',gA1,uip{:}),aup{:});
text(tx,ty,'A',tp{:});
axes('Parent',uicontainer('Parent',gA1,uip{:}),aup{:});
x=ctsbma(:,1);y=ctsbma(:,2);arg=[];arg.show=[0 1 1];
plotARP(x,y,arg);
xlabel('pre',afs{:}); 
ylabel('dur',afs{:},fw{:});
title('time correlations',tfs{:}); legend off;

set(gA1,'Heights', [25 -1]);


%% C1 : Before During Spatial
gC1 = uix.GridFlex('Parent',gTop,gfp{:});
axes('Parent',uicontainer('Parent',gC1,uip{:}),aup{:});
text(tx,ty,'C',tp{:});
axes('Parent',uicontainer('Parent',gC1,uip{:}),aup{:});
x=ctsbma(:,4);y=ctsbma(:,5);arg=[];arg.show=[0 1 1];
plotARP(x,y,arg);
xlabel('pre',afs{:}); 
ylabel('dur',afs{:},fw{:});
title('spatial correlations',tfs{:}); legend off;

set(gC1,'Heights', [25 -1]);

%% A2 : Before After Temporal
gA2 = uix.GridFlex('Parent',gTop,gfp{:});
axes('Parent',uicontainer('Parent',gA2,uip{:}),aup{:});
text(tx,ty,'',tp{:});
axes('Parent',uicontainer('Parent',gA2,uip{:}),aup{:});
x=ctsbma(:,1);y=ctsbma(:,3);arg=[];arg.show=[0 1 1];
plotARP(x,y,arg);
xlabel('pre',afs{:}); 
ylabel('post',afs{:},fw{:});

set(gA2,'Heights', [25 -1]);


%% C2 : Before After Spatial
gC2 = uix.GridFlex('Parent',gTop,gfp{:});
axes('Parent',uicontainer('Parent',gC2,uip{:}),aup{:});
text(tx,ty,'',tp{:});
axes('Parent',uicontainer('Parent',gC2,uip{:}),aup{:});
x=ctsbma(:,4);y=ctsbma(:,6);arg=[];arg.show=[0 1 1];
plotARP(x,y,arg);
xlabel('pre',afs{:}); 
ylabel('post',afs{:},fw{:});

set(gC2,'Heights', [25 -1]);

N = nshuffle;
%% HISTOGRAMs
gB = uix.GridFlex('Parent',gTop,gfp{:});
axes('Parent',uicontainer('Parent',gB,uip{:}),aup{:});
text(0,0.0,'B',tp{:});
gBB = uix.GridFlex('Parent',gB,gfp{:});
axes('Parent',uicontainer('Parent',gBB,uip{:}),aup{:});
bar([sum(ptsbma(:,1)<= N*0.01)/length(pairs), sum(N*0.99 <= ptsbma(:,1))/length(pairs);...
     sum(ptsbma(:,2)<= N*0.01)/length(pairs), sum(N*0.99 <= ptsbma(:,2))/length(pairs);...
     sum(ptsbma(:,3)<= N*0.01)/length(pairs), sum(N*0.99 <= ptsbma(:,3))/length(pairs);]*100,'stacked')
ylabel('% total',fss,10);
set(gca,'xticklabel',{'pre';'during';'post'},'YLim',[0,100],'xticklabelrotation',90);
title(sprintf('shuffling significance'));
legend({'(-)  corr';'(+) corr'}); axis square;
set(gB,'Heights', [25 -1]);

gD = uix.GridFlex('Parent',gTop,gfp{:});
axes('Parent',uicontainer('Parent',gD,uip{:}),aup{:});
text(0,0.0,'D',tp{:});
gDD = uix.GridFlex('Parent',gD,gfp{:});
axes('Parent',uicontainer('Parent',gDD,uip{:}),aup{:});
bar([sum(ptsbma(:,4)<= N*0.01)/length(pairs), sum(N*0.99 <= ptsbma(:,4))/length(pairs);...
     sum(ptsbma(:,5)<= N*0.01)/length(pairs), sum(N*0.99 <= ptsbma(:,5))/length(pairs);...
     sum(ptsbma(:,6)<= N*0.01)/length(pairs), sum(N*0.99 <= ptsbma(:,6))/length(pairs);]*100,'stacked')
ylabel('% total',fss,10);
set(gca,'xticklabel',{'pre';'during';'post'},'YLim',[0,100],'xticklabelrotation',90);
title(sprintf('shuffling significance'));
legend({'(-)  corr';'(+) corr'}); axis square;

set(gD,'Heights', [25 -1]);
%space

%% epilogue
set(gTop,'Widths', [-1 -1 -1.5]);%,'Widths', [-1 -1 -1]);

clear gTop; clear gA1;clear gA2;clear gC1;clear gC2; clear gB; clear gD;   

