function CellTab(fig, tab, m)
gui.Window = fig;
gui.m = m;
gui.m.sesh = 'before';
gui.m.ind = 1;
gui.m.cgi = [];
g = gui.m.groups{1};
gui.m.c = g(1);
for i = 1:length(gui.m.groups)
    g = gui.m.groups{i};
    for j = 1:length(g)
        c = g(j);
        gui.m.cgi(c.ind,:) = [i j] ;
    end
end
s = setSesh(gui.m.ind);
gui.m.t0 = 1;
gui.m.tn = length(gui.m.c.pt);

%TOP UI
CellTabBox = uix.HBoxFlex( 'Parent', tab, 'Spacing', 3);
controlPanel = uix.Panel('Parent', CellTabBox,'Title', '' );
gui.viewPanel = uix.Panel('Parent', CellTabBox,'Title', 'Ready','fontweight','bold',...
    'FontSize', 11, 'TitlePosition','centerTop','ForegroundColor','b');
set( CellTabBox, 'Widths', [200 -1] );

%% VIEW PANEL
gui.viewContainer = uicontainer('Parent', gui.viewPanel,'backgroundcolor','k');
%% CONTROL PANEL
controlBoxV = uix.VBox( 'Parent', controlPanel, 'Spacing', 3);
t = [];
%% CELL
cellsBoxV = uix.VBox( 'Parent', controlBoxV,'Padding', 3, 'Spacing', 0 ); t = [t -1];
uicontrol('Style','text','Parent', cellsBoxV,'HorizontalAlignment', 'left','String', 'Cells:');
gui.cellsListBox = uicontrol( 'Style', 'list','BackgroundColor', 'w','Parent', cellsBoxV, ...
    'String', num2str((1:length(gui.m.cgi))'),'Value', 1,'Callback', @onCellsListSelection);
set( cellsBoxV, 'Heights', [20 -1] );
%% MID
midBoxV = uix.VBox( 'Parent', controlBoxV,'Padding', 3, 'Spacing', 3 ); t = [t 150];
uicontrol('Parent', midBoxV,'Style','text','HorizontalAlignment', ...
    'left', 'String','Muscimol Bin:');
gui.midListBox = uicontrol( 'Style', 'list','BackgroundColor', 'w','String', {'before';'midall'}, ...
    'Parent', midBoxV, 'Value', 1,'Callback', @onMidListSelection, 'Enable', 'on');
set( midBoxV, 'Heights', [20 -1] );

%% t0
t0BoxV = uix.VBox( 'Parent', controlBoxV,'Padding', 3, 'Spacing', 3 , 'Visible', 'on'); t = [t 60];
gui.t0Label = uicontrol('Parent', t0BoxV,'Style','text',...
    'HorizontalAlignment', 'left','String', sprintf('t0 (s):%s',''));
t0BoxH = uix.HBox( 'Parent', t0BoxV,'Padding', 3, 'Spacing', 3 );
mnx = [1 floor(length(gui.m.c.pt)*0.02)]; step = 1; step = step/(mnx(2)-mnx(1)); %max-min
gui.t0Slider = uicontrol( 'Style', 'slider', 'Min', mnx(1),'Max', mnx(2), ...
    'Parent', t0BoxH,'Value', 1,'SliderStep',[step 10*step],'Callback', @ont0Slide);
gui.t0Edit = uicontrol( 'Parent', t0BoxH, 'Style', 'edit', 'String', 1,'Callback', @ont0Edit);
addlistener(gui.t0Slider,'ContinuousValueChange',@(hObject, event) ont0Sliding(hObject, event));
set( t0BoxH, 'Widths', [-4 -1] );

%% tn
tnBoxV = uix.VBox( 'Parent', controlBoxV,'Padding', 3, 'Spacing', 3 , 'Visible', 'on'); t = [t 60];
gui.tnLabel = uicontrol('Parent', tnBoxV,'Style','text',...
    'HorizontalAlignment', 'left','String', sprintf('tn (s):%s',''));
tnBoxH = uix.HBox( 'Parent', tnBoxV,'Padding', 3, 'Spacing', 3 );
mnx = [1 floor(length(gui.m.c.pt)*0.02)]; step = 1; step = step/(mnx(2)-mnx(1)); %max-min
gui.tnSlider = uicontrol( 'Style', 'slider', 'Min', mnx(1),'Max', mnx(2), ...
    'Parent', tnBoxH,'Value', mnx(2), 'SliderStep',[step 10*step],'Callback', @ontnSlide);
gui.tnEdit = uicontrol( 'Parent', tnBoxH, 'Style', 'edit', 'String', mnx(2),'Callback', @ontnEdit);
addlistener(gui.tnSlider,'ContinuousValueChange',@(hObject, event) ontnSliding(hObject, event));
set( tnBoxH, 'Widths', [-4 -1] );

%% FINAL HEIGHTS
set(controlBoxV, 'Heights', t );

%% CALLBACKS
% CELL TAB
    function onCellsListSelection( src, ~ )
        set(gui.viewPanel,'ForegroundColor',[0,0,1], 'Title', '');
        t = cellstr(get( src, 'String' ));
        gui.m.ind = str2double(t(get( src, 'Value' )));
        %gui.m.sesh = 'before'; %PUT BACK IN
        s = setSesh(gui.m.ind);
        set(gui.viewPanel,'Title',s);
        set(gui.midListBox,'String',{'before';'midall';'after'},'Value', 1);
        if ~gui.m.g(1).after.exists
            set(gui.midListBox,'String', {'before';'midall'});
        end
        if gui.m.tn > length(gui.m.c.pt)
            gui.m.tn = length(gui.m.c.pt);
        end
        gui.m.t0 = 1;
        gui.m.tn = length(gui.m.c.pt);
        mnx = [1 floor(length(gui.m.c.pt)*0.02)]; step = 1; step = step/(mnx(2)-mnx(1)); %max-min
        set(gui.tnSlider,'Min', mnx(1),'Max', mnx(2),'Parent', tnBoxH,'Value',  floor(length(gui.m.c.pt)*0.02) );
        run()
    end
% MID
    function onMidListSelection( src, ~ )
        t = get( src, 'String');
        gui.m.sesh = t{get(src,'Value')};
        s = setSesh(gui.m.ind);
        gui.m.t0 = 1;
        gui.m.tn = length(gui.m.c.pt);
        mnx = [1 floor(length(gui.m.c.pt)*0.02)]; step = 1; step = step/(mnx(2)-mnx(1)); %max-min
        set(gui.tnSlider,'Min', mnx(1),'Max', mnx(2),'Parent', tnBoxH,'Value', floor(length(gui.m.c.pt)*0.02) );
        set(gui.viewPanel,'Title',s);
        run()
    end
% t0
    function ont0Slide( src, ~ )
        t = round2r(get(src, 'Value'), 1);
        gui.m.t0 = t/0.02; %to timesteps
        set(gui.t0Edit, 'String', t);
        run();
    end
    function ont0Edit(src,~)
        t = round2r(str2double(get(src,'String')), 1);
        gui.m.t0 = t/0.02; %to timesteps;
        set(gui.t0Slider, 'Value', t);
        run();
    end
    function ont0Sliding(hObject,~)
        t = round2r(get(hObject,'Value'),1);
        set(gui.t0Edit, 'String', t);
    end
% t0
    function ontnSlide( src, ~ )
        t = round2r(get(src, 'Value'), 1);
        gui.m.tn = t/0.02; %to timesteps
        set(gui.tnEdit, 'String', t);
        run();
    end
    function ontnEdit(src,~)
        t = round2r(str2double(get(src,'String')), 1);
        gui.m.tn = t/0.02; %to timesteps;
        set(gui.tnSlider, 'Value', t);
        run();
    end
    function ontnSliding(hObject,~)
        t = round2r(get(hObject,'Value'),1);
        set(gui.tnEdit, 'String', t);
    end



%% UTIL
    function s = setSesh(ind)
        gui.m.g = gui.m.groups{gui.m.cgi(ind,1)};
        i = gui.m.cgi(gui.m.ind,2);
        if strcmp(gui.m.sesh,'before')
            gui.m.c = gui.m.g(i).before; s = 'before';
        elseif strcmp(gui.m.sesh,'midall')
            gui.m.c = gui.m.g(i).midall; s = 'midall';
        elseif strcmp(gui.m.sesh,'after')
            gui.m.c = gui.m.g(i).after; s = 'after';
        else
            gui.m.c = gui.m.g(i).middle{gui.m.sesh}; s = ['mid' num2str(gui.m.sesh)];
        end
        s = sprintf('Cell %d %s',gui.m.ind, s);
    end
    

function f = gaussian2d(N,sigma)
  % N is grid size, sigma speaks for itself
 [x y]=meshgrid(round(-N/2):round(N/2), round(-N/2):round(N/2));
 f=exp(-x.^2/(2*sigma^2)-y.^2/(2*sigma^2));
 f=f./sum(f(:));
end

%% RUN %%
    function run()
        delete(findobj(tab,'type','axes'));
        gui.m.parent = gui.viewContainer;
        if gui.m.t0 > gui.m.tn
            gui.m.tn = tui.m.t0 + 1;
        end
        %c.ind     
        c = gui.m.c; dt = round(median(diff(c.pt)),3); %to ms
        %a = patchTrajectoryLinear(c.pt,c.px,c.py,dt,dt*1.9); %NO NEED NOW?
        %c.pt = a.t; c.px = a.x; c.py = a.y;
        pw = gui.m.t0:gui.m.tn; %pw = 1:length(c.pt); %DO PW HERE
        c.si = discretize(c.st, [-Inf, mean([c.pt(2:end)'; c.pt(1:end-1)']), +Inf]);
        sw = c.si(c.si>gui.m.t0 & c.si< gui.m.tn);
        
        tfs = 16;
        rows = 3;
        cols = 3;
        axi = []; ms  = 7; lw = 2; 
        
        %PLOT 1        
        ax = subplot(rows,cols,1,'Parent',gui.m.parent); axi(end+1) = ax;
        %XY
        %set(ax,'color','w');
        %plot(ax,c.px(pw),c.py(pw),'b','linewidth',1);hold(ax,'on');
        plot(ax,c.px(sw),c.py(sw),'w.','markersize',3);%,'linewidth',ss);
        title(ax,sprintf('spikes=%d',length(sw)),'color','m','fontsize',tfs);
        
        %PLOT 2
        ax = subplot(rows,cols,8,'Parent',gui.m.parent); axi(end+1) = ax;
         %RM
        nbins = 50;
        mx = max(c.px(pw));
        my = max(c.py(pw));
        pxi = discretize(c.px(pw), 0:mx/nbins:mx);
        pyi = discretize(c.py(pw), 0:my/nbins:my);
        %t = diff(c.pt(pw)); t = [median(t); t];
        rmt = accumarray([pyi pxi], 1, [nbins nbins]);
        %rmt = accumarray([pyi pxi], t, [nbins nbins]);
        rmt(rmt<1) = 1e10;
        %rmt(rmt<min(t)) = 1e10;
        %rmt = accumarray([pxi' pyi'], 1, [nbins nbins]); %SORT OUT WITH FUNCTION
        sxi = discretize(c.px(sw), 0:mx/nbins:mx);
        syi = discretize(c.py(sw), 0:my/nbins:my);
        rms = accumarray([syi sxi], 1, [nbins nbins]);
        %rms = accumarray([sxi' syi'], 1, [nbins nbins]);
        rm = rms./(rmt*dt);
        %rm = rms./rmt;
        [m, i] = max(rm(:));
        ind2sub(size(rm),i);
        %t = rms-rmt; t(t<0)=0; rmt = rmt + t; %add in extra timestep ONLY when spikes occured faster than timestep
        rm = (rms ./ (rmt*dt));
        rm(isnan(rm)) = 0; %DO THIS?
        %maxfiring rate
        %rm = imgaussfilt(rm,1);
        %rm = rm/max(rm(:)); %normalize %REMOVE??        
        %f = gaussian2d(nbins,2);
        %imagesc(ax, conv2(rm,f,'same'));
        imagesc(ax, imgaussfilt(c.rm,1)); %rm
        title(ax,sprintf('max %.1fHz',m),'color','m','fontsize', tfs);
        
        %PLOT 3
        ax = subplot(rows,cols,3,'Parent',gui.m.parent); axi(end+1) = ax;
        % AC
        ac = xcorr2( rm);
        ac = ac/ max(ac(:)); %normalized ac;
        acg = imgaussfilt(ac, 2,'FilterDomain','spatial'); %imshow
        imagesc(ax,acg); hold(ax,'on'); %%%%
        [k,l] = find(imregionalmax(acg));
        dist = pdist2([l k],[nbins, nbins] ); [dist,ind]=sort(dist);
        if length(l) >= 7
            l=l(ind(2:7)); k=k(ind(2:7)); 
        end
        %plot(ax,l,k,'m+','markersize',s,'linewidth',ss);
        acg = imgaussfilt(ac, 3,'FilterDomain','spatial'); %imshow
        [k,l] = find(imregionalmax(acg));
        dist = pdist2([l k],[nbins, nbins] ); [dist,ind]=sort(dist);
        if length(l) >= 7
            l=l(ind(2:7)); k=k(ind(2:7)); 
        end
        plot(ax,l,k,'yx','markersize',ms,'linewidth',lw);
        %plot(ax,[nbins l(2:3)' nbins],[nbins k(2:3)' nbins],'m'); triangle
        if length(l) >= 2
            [~,ceny]=max(max(ac));[~,cenx]=max(max(ac'));
            %viscircles(ax,[cenx ceny], dist(2)/2,'color','m'); %ADD BACK
            %viscircles(ax,[cenx ceny], dist(2)*3/2,'color','m');
        end
        %rms(rms>0)=1;
        title(ax,sprintf('gs %.2f gs2 %.2f gs3 %.2f',gridscore2(ac,2),...
            gridscore2(xcorr2(imgaussfilt(rm,2)),2),...
            gridscore2(ac,3)),'color','m','fontsize',tfs); %rms
        % ax = subplot(rows,cols,3,'Parent',gui.m.parent)  %PLOT 1 TITLE
         %title(ax,sprintf('gs %.2f \ngs2 %.2f \ngs3 %.2f',gridscore2(ac,2),...
             %gridscore2(xcorr2(imgaussfilt(rms,1)),2),...
             %c.gridscore),'color','m','fontsize',18);
        % AC2 SMOOTHED
        %ax = subplot(rows,cols,7,'Parent',gui.m.parent) ; axi(end+1) = ax;
        %[cs, rs] = imfindcircles(ac>std(ac(:)),[10 40]);
        %imagesc(ax,acg); %hold on;
        
	     % $$ STILL ADDING TIME TO BINS
        % $$ STILL ADDING TIME TO BINS
        % 36 12.0
        
        % PLOT 4
        ax = subplot(rows,cols,4,'Parent',gui.m.parent);axi(end+1) = ax;
        %MODULE
        acg = imgaussfilt(ac, 3,'FilterDomain','spatial');
        %acg = imresize(acg, [70 120]);
        imagesc(ax,acg); hold(ax,'on'); %%%%
        cent = size(acg)/2;
        plot(ax,cent(2),cent(1),'wx','markersize',10);
        module = Find_Module(acg);
        plot(ax,module.hex_peaks(:,1),module.hex_peaks(:,2),'ro');
        %plot(ax,module.hex_peaks2(:,1),module.hex_peaks2(:,2),'y*'); 
        plot(ax,module.x,module.y,'r'); 
        title(ax,sprintf(' %s','Module'),'color','m','fontsize',tfs);
        
        % PLOT 5
         ax = subplot(rows,cols,5,'Parent',gui.m.parent); axi(end+1) = ax;
         % RAYLEIGH
         d = gui.m.c;
         rd = histcounts(rad2deg(c.hd(c.si))+180,0:3:360)./histcounts(rad2deg(c.hd)+180,0:3:360);
         rd = smooth(rd,15);
         plot(ax, 3:3:360,rd,'g','linewidth',1.8);
         rs = d.rayleigh_score;
         title(ax,sprintf('HD Rate rayleigh=%.2f',rs),'color','m','fontsize',tfs);
        
         % PLOT 6
         ax = subplot(rows,cols,6,'Parent',gui.m.parent); axi(end+1) = ax;
         % TIME DIFF
         plot(ax, c.pt); 
         hold(ax,'on'); plot(ax,c.si,c.st,'r.');
         df = diff(c.pt);
         hold(ax,'off'); hist(ax,df,10);
         %title(ax,sprintf('Time diff mx%.2f m%.3f s%.2f', max(df), mean(df), std(df) ),'color','m','fontsize',tfs-2);
         title(ax,sprintf('%s','Time diff (s)'),'color','m','fontsize',tfs);
         
         % PLOT 7
         ax = subplot(rows,cols,7,'Parent',gui.m.parent); axi(end+1) = ax;
         % TRAIN
         t = zeros(c.si(end),1);
         t(c.si)=1;
         bin_s = 2;
         t = histcounts(c.st,0:bin_s:c.st(end));
         plot(ax,t,'y');
         title(ax,sprintf('Train bin=%.2fs', bin_s),'color','m','fontsize',tfs);
         
         
         
         % FOR PAPER - REMOVE
         ax = subplot(rows,cols,2,'Parent',gui.m.parent); axi(end+1) = ax; cla(ax);
         % AC Clean
         acg = imgaussfilt(ac, 3,'FilterDomain','spatial');imagesc(ax,acg);
         title(ax,sprintf('gridscore=%.1f',round(gridscore2(ac,3),1)),'color','m','fontsize',tfs); %rms
         
         ax = subplot(rows,cols,3,'Parent',gui.m.parent); axi(end+1) = ax; cla(ax);
         % RAYLEIGH
         d = gui.m.c;
         rd = histcounts(rad2deg(c.hd(c.si))+180,0:3:360)./histcounts(rad2deg(c.hd)+180,0:3:360);
         rd = smooth(rd,45);
         plot(ax, 3:3:360,rd,'g','linewidth',3);
         rs = d.rayleigh_score;
         title(ax,sprintf('Rayleigh=%.1f',round(rs,1)),'color','m','fontsize',tfs); 
         xlabel(ax,'angle','color','m','fontsize',tfs-2)
        
        %ALLLL
        for i = 1:length(axi)
            colormap(axi(i), 'jet');
            set(axi(i),'YDir','normal');
            axis(axi(i),'tight');
            set(axi(i),'Ycolor','m');
            set(axi(i),'Xcolor','m');
            set(axi(i),'Color','k');
            axis(axi(i),'square');
        end
        
    end %end run()


end








        
        
        %center
        %ax = subplot(rows,cols,9,'Parent',gui.m.parent); axi(end+1) = ax;
%         plot(ax,c.px,c.py,'linewidth',1);hold(ax,'on');
%         plot(ax,c.px(c.si),c.py(c.si),'w.','markersize',8);%,'linewidth',ss);
%         set(ax,'Color','k');
        
        
%         % RM PREV
%         ax = subplot(rows,cols,3,'Parent',gui.m.parent); axi(end+1) = ax;
%         imagesc(ax,c.rm);hold(ax,'on');
%         
%         % AC PREV
%         ac = c.ac;
%         ax = subplot(rows,cols,6,'Parent',gui.m.parent); axi(end+1) = ax;
%         imagesc(ax,ac);hold(ax,'on');  %%%%
%         ac(isnan(ac))=0;
%         [k,l] = find(imregionalmax(ac));
%         plot(ax,l,k,'m+','markersize',s,'linewidth',ss);
%         [zmax,imax,zmin,imin]= Extrema2(ac);
%         [i,j]=ind2sub(size(ac),imax); %NOAM
%         %plot(ax,j,i,'mx','markersize',s,'linewidth',ss);
%         title(ax,sprintf('gs %.2f',gridscore(ac,-1)),'color','m','fontsize',18);
        
        

  
%         %AC Cross_Correlation
%         ac = Cross_Correlation(rm, rm);
%         ax = subplot(rows,cols,5,'Parent',gui.m.parent); axi(end+1) = ax;
%         ac(isnan(ac)) = 0;
%         %ac(nbins,nbins) = max(ac(:))+0.01; %for gridscore function
%         aco = ac;
%         ac = conv2(ac,gaussian2d(length(ac),2),'same');
%         imagesc(ax,aco);hold(ax,'on');  %%%%
%         [k,l] = find(imregionalmax(ac));
%         %put all extrema points in dist
%          ms = min(size(ac))/2; mss = ms/4; mss = round((-mss:1:mss) + ms);t = ac(mss,mss);
%         [cenx ceny]=ind2sub(size(ac),find(ac==max(t(:))));
%         dist = pdist2([l k],[cenx, ceny] ); [dist,ind]=sort(dist);
%         if length(l) >= 7
%             l=l(ind(2:7)); k=k(ind(2:7)); 
%         end
%         plot(ax, l,k,'m+','markersize',s,'linewidth',ss);
%         viscircles(ax,[cenx ceny], dist(2)/2,'color','m');
%         viscircles(ax,[cenx ceny], dist(7),'color','m');
%         title(ax,sprintf('gs %.2f',gridscore2(aco,-1)),'color','m','fontsize',18);
%         %AC CC SMOOTHED 
%         ax = subplot(rows,cols,8,'Parent',gui.m.parent); axi(end+1) = ax;
%         %ac = imgaussfilt(aco, 2,'FilterDomain','spatial'); %imshow same
%         imagesc(ax,ac);hold(ax,'on');


        

               
        %imagesc(ax,imgaussfilt(rm, 1)); %imshow
        % CORRELATION
        %imagesc(corrcoef(rm));
        %imagesc(xcorr2(rm,rm));
        %a(a > 0.1) = 0.1; 
        %imagesc(a);



