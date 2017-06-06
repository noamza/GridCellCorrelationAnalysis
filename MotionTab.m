function MotionTab(fig, tab, m)
    gui.Window = fig;
    gui.m = m;
    gui.m.speed = 200;
    gui.m.sesh = 'before';
    gui.m.cells = ones(length(gui.m.g),1);
    gui.m.pause = false;
    gui.m.stop = false;
    gui.m.go = 0;
    gui.m.timestep = 1;
    gui.m.win = 500;
    gui.m.showTrain = false;
    %gui.m.movie = VideoWriter('gridmotion.avi'); open(gui.m.movie);
    gui.m.record = false;
    updateM();
    
    
    %TOP UI
    AniTabBox = uix.HBoxFlex( 'Parent', tab, 'Spacing', 3);
    controlPanel = uix.Panel('Parent', AniTabBox,'Title', '' );
    gui.viewPanel = uix.Panel('Parent', AniTabBox,'Title', 'Ready','fontweight','bold',...
        'FontSize', 11, 'TitlePosition','centerTop','ForegroundColor','b');
    set( AniTabBox, 'Widths', [200 -1] );
    
    %% VIEW PANEL
    gui.viewContainer = uicontainer('Parent', gui.viewPanel,'backgroundcolor','k');
    
    %% CONTROL PANEL
    controlBoxV = uix.VBox( 'Parent', controlPanel, 'Spacing', 3);
    t = [];
    %GROUP
    groupBoxV = uix.VBox( 'Parent', controlBoxV,'Padding', 3, 'Spacing', 0 ); t = [t -1];
    uicontrol('Style','text','Parent', groupBoxV,'HorizontalAlignment', 'left','String', 'Group:');
    gui.groupListBox = uicontrol( 'Style', 'list','BackgroundColor', 'w','Parent', groupBoxV, ...
        'String', num2str((1:length(gui.m.groups))'),'Value', 1,'Callback', @onGroupListSelection);
    set( groupBoxV, 'Heights', [20 -1] );
    %MID
    midBoxV = uix.VBox( 'Parent', controlBoxV,'Padding', 3, 'Spacing', 3 ); t = [t 150];
    uicontrol('Parent', midBoxV,'Style','text','HorizontalAlignment', ...
        'left', 'String','Muscimol Bin:');
    gui.midListBox = uicontrol( 'Style', 'list','BackgroundColor', 'w', ...
        'Parent', midBoxV,'String', {'before';'midall';'after'}, ...
        'Value', 1,'Callback', @onMidListSelection, 'Enable', 'on');
    set( midBoxV, 'Heights', [20 -1] );
    %CHECK BOXES
    gui.checkBoxV = uix.VBox('Parent', controlBoxV); t = [t 250];
    uicontrol('Style','text','Parent', gui.checkBoxV, 'String', 'Show Cells:',...
        'HorizontalAlignment', 'left');
    gui.cbh = zeros(length(gui.m.g),1);
    for i = 1:length(gui.m.g)
        c = getSesh(i);
        s = sprintf('i%d tet%d grid%.1f', gui.m.g(i).ind, gui.m.g(i).tet, c.gridscore);
        gui.cbh(i) = uicontrol('Style','checkbox','String', s,'Value',1,...
                'Parent', gui.checkBoxV, 'Callback',{@checkBoxCallback,i},...
                'ForegroundColor',gui.m.colors(i,:), 'FontSize', 10, 'fontweight', 'bold');
    end
     %GO
    goBoxH = uix.HBox( 'Parent', controlBoxV, 'Spacing', 3); t = [t 30];
    uicontrol('Parent', goBoxH,'Style','text', 'String', 'Go to (min): ',...
        'HorizontalAlignment', 'left');
    gui.goEdit = uicontrol( 'Parent', goBoxH, 'Style', 'edit', 'String', gui.m.go, ...
        'Callback', @onGoEdit, 'FontSize', 9);
    set( goBoxH, 'Widths', [125 -1] );
    %TIME
    timeBoxH = uix.VBox( 'Parent', controlBoxV,'Padding', 3, 'Spacing', 3 , 'Visible', 'on'); t = [t 60];      
    uicontrol('Parent', timeBoxH,'Style','text',...
        'HorizontalAlignment', 'left','String', 'Time Step:');
    timeBoxH = uix.HBox( 'Parent', timeBoxH,'Padding', 3, 'Spacing', 3 );
    mnx = [1 1000]; step = 1; step = step/(mnx(2)-mnx(1)); %max-min 
    gui.timeSlider = uicontrol( 'Style', 'slider', 'Min', mnx(1),'Max', mnx(2), ...
        'Parent', timeBoxH,'Value', gui.m.timestep,'SliderStep',[step 10*step],'Callback', @onTimeSlide);        
    gui.timeEdit = uicontrol( 'Parent', timeBoxH, 'Style', 'edit', 'String', gui.m.timestep,...
        'Callback', @onTimeEdit);
    addlistener(gui.timeSlider,'ContinuousValueChange',@(hObject, event) onTimeSliding(hObject, event));
    set( timeBoxH, 'Widths', [-4 -1] );
    %SPEED
    speedBoxH = uix.VBox( 'Parent', controlBoxV,'Padding', 3, 'Spacing', 3 , 'Visible', 'on'); t = [t 60];      
    uicontrol('Parent', speedBoxH,'Style','text','String', 'Speed (lightyears):',...
        'HorizontalAlignment', 'left');
    speedBoxH = uix.HBox( 'Parent', speedBoxH,'Padding', 3, 'Spacing', 3 );
    mnx = [1 2000]; step = 10; step = step/(mnx(2)-mnx(1)); %max-min 
    gui.speedSlider = uicontrol( 'Style', 'slider', 'Min', mnx(1),'Max', mnx(2), ...
        'Parent', speedBoxH,'Value', gui.m.speed,'SliderStep',[step 3*step],'Callback', @onSpeedSlide);        
    gui.speedEdit = uicontrol( 'Parent', speedBoxH, 'Style', 'edit', 'String', gui.m.speed,...
        'Callback', @onSpeedEdit);
    addlistener(gui.speedSlider,'ContinuousValueChange',@(hObject, event) onSpeedSliding(hObject, event));
    set( speedBoxH, 'Widths', [-4 -1] );
    %WINDOW
    winBoxV = uix.VBox( 'Parent', controlBoxV,'Padding', 3, 'Spacing', 3 , 'Visible', 'on'); t = [t 60];
    winBoxH1 = uix.HBox( 'Parent', winBoxV,'Padding', 3, 'Spacing', 3 );
    gui.showTrainCh = uicontrol('Style','checkbox','Value',gui.m.showTrain,'String', 'Train Length(s):',...
          'Parent', winBoxH1, 'Callback',@OnshowTrainCh); %'FontSize', 10, 'fontweight', 'bold'
    %uicontrol('Parent', winBoxH1,'Style','text','String', 'Train Window(s):', 'HorizontalAlignment', 'left');
    winBoxH = uix.HBox( 'Parent', winBoxV,'Padding', 3, 'Spacing', 3 );
    mnx = [50 20000]; step = 50; step = step/(mnx(2)-mnx(1)); %max-min 
    gui.winSlider = uicontrol( 'Style', 'slider', 'Min', mnx(1),'Max', mnx(2), ...
        'Parent', winBoxH,'Value', gui.m.win,'SliderStep',[step 3*step],'Callback', @onWinSlide);        
    gui.winEdit = uicontrol( 'Parent', winBoxH, 'Style', 'edit', 'String', gui.m.win,...
        'Callback', @onWinEdit);
    addlistener(gui.winSlider,'ContinuousValueChange',@(hObject, event) onWinSliding(hObject, event));
    set( winBoxH, 'Widths', [-4 -1] );
    %BUTTONS
    buttonBoxH = uix.HBox( 'Parent', controlBoxV, 'Spacing', 3); t = [t 50];
    gui.RunButton = uicontrol( 'Style', 'PushButton','Parent', buttonBoxH,...
        'String', 'Start','Callback', @run );
    gui.PauseButton = uicontrol( 'Style', 'PushButton','Parent', buttonBoxH,...
        'String', 'Pause','Callback', @onPause,'Enable','On');
    gui.StopButton = uicontrol( 'Style', 'PushButton','Parent', buttonBoxH,...
        'String', 'Stop','Callback', @onStop,'Enable','On');
    gui.RecordButton = uicontrol( 'Style', 'PushButton','Parent', buttonBoxH,...
        'String', 'Record','Callback', @onRecord,'Enable','On','backgroundColor','w');
    %STATUS
    %gui.status = uicontrol('Parent', controlBoxV,'Style','text', 'Tag', 'status','String', 'listo',...
    %'HorizontalAlignment', 'left', 'FontSize',9,'foregroundcolor','b'); t = [t 30];
    %FINAL HEIGHTS
    set(controlBoxV, 'Heights', t );
    
    function makeColors()
        disp('colors');
        gui.m.colors = jet(length(gui.m.g)); %jet hsv
    end
    
    %% CALLBACKS
    %CHECKBOXES
    function checkBoxCallback(hObject,eventData,id)
        gui.m.cells(id) = get(hObject,'Value'); [id get(hObject,'Value')] 
        gui.m.cells
    end
    %GROUP
    function onGroupListSelection( src, ~ )
        %set(gui.status, 'ForegroundColor',[0,0,1]);
        onStop(-1,-1);
        set(gui.viewPanel,'ForegroundColor',[0,0,1], 'Title', '');
        t = cellstr(get( src, 'String' ));
        gui.m.gid = str2double(t(get( src, 'Value' )));
        gui.m.g = gui.m.groups{gui.m.gid};
        updateM();
        l = length(gui.m.g);
        disp(['g: ', num2str(gui.m.gid),' ', num2str(l)]);
        set(gui.viewPanel,'Title', sprintf('group size: %d',l));
        %set(gui.status, 'String', sprintf('group size: %d',l));cellstr((num2str(1:length(gui.m.g(1).middle),'%d'))');
        set(gui.midListBox,'String',['before';'midall';'after'],'Value', 1);
        %reset checkboxes
        delete(findobj(gui.checkBoxV.Children,'Style','checkbox'));
        gui.cbh = zeros(length(gui.m.g),1);
        gui.m.cells = ones(length(gui.m.g),1);
        for i = 1:length(gui.m.g)
            c = getSesh(i);
            s = sprintf('i%d tet%d grid%.1f', gui.m.g(i).ind, gui.m.g(i).tet, c.gridscore);
            gui.cbh(i) = uicontrol('Style','checkbox','String', s,'Value',1,...
                'Parent', gui.checkBoxV, 'Callback',{@checkBoxCallback,i},...
                'ForegroundColor',gui.m.colors(i,:), 'FontSize', 10, 'fontweight', 'bold');
        end
        [~,s] = getSesh(1);
        set(gui.viewPanel,'Title', s);
        gui.m.go = 1; %$$$$$$$$$$$$$$$
        run();  %REDO THIS
    end
    % MID
    function onMidListSelection( src, ~ )
        %gui.m.mid = get( src, 'Value' )-1;
        gui.m.stop = true;
        t = get( src, 'String');
        gui.m.sesh = str2double(t{get( src, 'Value' )});
        if isnan(gui.m.sesh)
            gui.m.sesh = t{get(src,'Value')};
        end
        updateM();
        %checkboxes
        for i = 1:length(gui.m.g)
            c = getSesh(i);
            s = sprintf('i%d tet%d grid%.1f', gui.m.g(i).ind, gui.m.g(i).tet, c.gridscore);
            set(gui.cbh(i),'String', s);
        end
        [~,s] = getSesh(1);
        set(gui.viewPanel,'Title', s);
        run(); %REDO THIS
        %}
    end
    %STOP
    function onStop(~, ~)
        gui.m.stop = true;
        gui.m.go = 1;
    end
    %PAUSE
    function onPause(~, ~)
        %delete(findobj(gui.Window.Children.Children(gui.Window.Children.Selection),...
        %'type','axes'));
        if gui.m.pause
            set(gui.PauseButton, 'String', 'Pause')
        else
            set(gui.PauseButton, 'String', 'Resume')
        end
        drawnow
        gui.m.pause = ~gui.m.pause;
    end
    %Go
    function onGoEdit(src,~)
        c = getSesh(1);
        gui.m.go = round(c.pt(1) + str2double(get(src,'String'))*60/.02 );
        gui.m.go;
    end
   function onRecord(~, ~)
        if gui.m.record %stop
            %figure;
            %movie(gui.m.movie);
            close(gui.m.movie);
            set(gui.RecordButton, 'String', 'Record','backgroundColor','w')
        else %start
            %gui.m.movie = VideoWriter(sprintf('Grid_Motion_%s.mp4',datestr(now,'yyyy.mm.dd_HH.MM.SS'))); %','
            gui.m.movie = VideoWriter(sprintf('Grid_Motion_%s.avi',datestr(now,'yyyy.mm.dd_HH.MM.SS'))); 
            %gui.m.movie.Quality = 100;
            open(gui.m.movie);
            set(gui.RecordButton, 'String', 'Record','backgroundColor','r');
        end
        gui.m.record = ~gui.m.record ;
    end
    %TIME
    function onTimeSlide( src, ~ )
        t = round2r(get(src, 'Value'), 1);
        gui.m.timestep = t;
        set(gui.timeEdit, 'String', t);
        %run();
    end
    function onTimeEdit( src, ~ )
        t = round2r(str2double(get(src, 'String')),1);
        gui.m.timestep = t;
        set(gui.timeSlider, 'Value', t);
    end
    function onTimeSliding(hObject,event)
        t = round2r(get(hObject,'Value'),1);
        set(gui.timeEdit, 'String', t);
    end 
    %SPEED
    function onSpeedSlide( src, ~ )
        t = round2r(get(src, 'Value'), 1);
        gui.m.speed = t;
        set(gui.speedEdit, 'String', t);
        %run();
    end
    function onSpeedEdit( src, ~ )
        t = round2r(str2double(get(src, 'String')),1);
        gui.m.speed = t;
        set(gui.speedSlider, 'Value', t);
    end
    function onSpeedSliding(hObject,event)
        t = round2r(get(hObject,'Value'),1);
        set(gui.speedEdit, 'String', t);
    end 
%WINDOW
function OnshowTrainCh(src, ~ )
    gui.m.showTrain = get(src, 'Value');
end
function onWinSlide( src, ~ )
        t = round2r(get(src, 'Value'), 100);
        gui.m.win = t;
        set(gui.winEdit, 'String', t);
        gui.m.lx = linspace(0,100,gui.m.win);
        gui.m.ly = repmat((length(gui.m.g):-1:1)',1,gui.m.win)*2 + max(gui.m.c.py) + 2;
        %run();
    end
    function onWinEdit( src, ~ )
        t = round2r(str2double(get(src, 'String')),100);
        gui.m.win = t;
        gui.m.lx = linspace(0,100,gui.m.win);
        gui.m.ly = repmat((length(gui.m.g):-1:1)',1,gui.m.win)*2 + max(gui.m.c.py) + 2;
        set(gui.winSlider, 'Value', t);
    end
    function onWinSliding(hObject,event)
        t = round2r(get(hObject,'Value'),100);
        set(gui.winEdit, 'String', t);
    end 
    
    
    function updateM()
        gui.m.g = gui.m.groups{gui.m.gid};
        gui.m.pause = false;
        gui.m.stop = false; 
        gui.m.c = getSesh(1);
        makeColors();
        train();
        gui.m.lx = linspace(0,100,gui.m.win);
        gui.m.ly = repmat((length(gui.m.g):-1:1)',1,gui.m.win)*2 + max(gui.m.c.py) + 2;
    end
    
    %% RUN %%
    function run(~,~)
        disp('RUN');
        delete(findobj(tab,'type','axes'));
        [~,s] = getSesh(1);
        set(gui.viewPanel,'Title', s);
        gui.m.parent = gui.viewContainer;
        gui.m.go = gui.m.go + 1; %so not 0;
        %train();
        z = -1;
        runPrivate(gui.m, z);
    end %end run()
    
    function runPrivate(m, z)
        %set(gui.status, 'String', sprintf('group(%d) session(%s)',gui.m.gid, t));
        if z == -1
            z = length(gui.m.c.pt);
        end
        ax = axes('Parent',gui.m.parent);
        %figure; ax = gca; %DELETE
        xlim(ax,[0 120]); ylim(ax,[0 120])
        %DRAW LOOP
        hold(ax,'off');
        axis(ax,'off');
        i = 1;
        %tic
        %START WHILE LOOP
        c = gui.m.c;
        mf = 1;
        xoff = 5; yoff = 5;
        c.px = c.px + xoff;c.py = c.py + yoff;
        maxx = max(c.px);
        maxy = max(c.py); 
        while gui.m.go <= z
            if gui.m.stop
                gui.m.stop = false;
                cla(ax);
                break
            end
            while(gui.m.pause)
                pause(1/5);
            end
            if gui.m.go > length(c.pt)
                gui.m.go = max(c.pt(1),length(c.pt) - round(1/0.02)); %set it back 1 min from end
            end
            i = gui.m.go;
            plot(ax,c.px(1:i),c.py(1:i),'Color',[0.3 0.3 0.3],'linewidth', 3);
            hold(ax, 'on');
            xlim(ax, [0 maxx+20]); ylim(ax, [0  maxy+20]);
            plot(ax,[xoff,xoff,maxx,maxx,xoff],[yoff,maxy+1,maxy+1,yoff,yoff],'w','linewidth', 2);
            axis(ax,'off');
            %plot cells
            tx = c.px(1:i); ty = c.py(1:i);
            %l = randperm(length(gui.m.g));
            for j = 1:length(gui.m.g)
                %if train(j,i) && m.cells(j)
                if gui.m.cells(j)
                    plot(ax, double(tx(gui.m.train(j,1:i))), double(ty(gui.m.train(j,1:i)))+0.3*j,...
                        '.','Color', gui.m.colors(j,:),'MarkerSize',15); %'markerfacecolor',gui.m.colors(j,:)
                    %plot(ax,lx(gui.m.train(j,1:i)),ly(j,gui.m.train(j,1:i)),'.','Color', colors(j,:))
                    wine = max(1,i-gui.m.win+1);
                    %plot(ax,lx(gui.m.train(j,wine:i)),ly(j,gui.m.train(j,wine:i)),'.','Color', colors(j,:));
                    if(gui.m.showTrain)
                        plot(ax,gui.m.lx(1:i-wine+1)+xoff,gui.m.ly(j,1:i-wine+1)+gui.m.train(j,wine:i)*1.9+yoff,...
                            'Color', gui.m.colors(j,:));
                        plot(ax,gui.m.lx(1:i-wine+1)+xoff,gui.m.ly(j,1:i-wine+1)+yoff,'Color', [0.2 0.2 0.2]);
                    end
                end
            end
            text(ax, double(maxx+2), double(maxy+5+yoff),sprintf('%.2f',c.pt(double(i))-c.pt(1)),...
                'fontsize',14,'color','w');
            hold(ax, 'off');
            drawnow();
            if(gui.m.record)
                %gui.m.movie(mf) = getframe(ax);
                frame = getframe(ax);
                writeVideo(gui.m.movie,frame);
                mf = mf+1;
            end
            %delete(ax);
            pause(1/gui.m.speed);
            %i = i + 1;
            gui.m.go = gui.m.go + gui.m.timestep;
        end
        %toc
    end%runPriv()
    
    function train
        train = [];%zeros(length(m.g), length(c.pt));
        c = [];
        for i = 1:length(gui.m.g)
            c = getSesh(i);
            if strcmp(gui.m.sesh,'before')
                edges = [-Inf, mean([c.pt(2:end)'; c.pt(1:end-1)']), +Inf];
            elseif strcmp(gui.m.sesh,'midall')
                edges = [-Inf, mean([c.pt(2:end); c.pt(1:end-1)]), +Inf];
            elseif strcmp(gui.m.sesh,'after')
                edges = [-Inf, mean([c.pt(2:end)'; c.pt(1:end-1)']), +Inf];
            else
                c = gui.m.g(i).middle{gui.m.sesh};
                edges = [-Inf, mean([c.pt(2:end); c.pt(1:end-1)]), +Inf];
            end
            %if isnumeric(gui.m.sesh)
            if isnan(edges(2))
                edges(2) = 1.1;
            end
            I = discretize(c.st, edges);
            train(i,:) = zeros(1, length(c.pt));
            train(i,I) = 1;
        end
        gui.m.train = logical(train);
    end
    
    function [c s] = getSesh(i)
        if strcmp(gui.m.sesh,'before')
            c = gui.m.g(i).before; s = 'before';
    elseif strcmp(gui.m.sesh,'midall')
            c = gui.m.g(i).midall; s = 'midall';
        elseif strcmp(gui.m.sesh,'after')
            c = gui.m.g(i).after; s = 'after';
        else
            c = gui.m.g(i).middle{gui.m.sesh}; s = ['mid' num2str(gui.m.sesh)];
        end
        s = sprintf('group %d session %s',gui.m.gid, s);
    end
end




function test(c, z)
    i = 0;
    figure;
    hold off
    tic
    while i < z
        i = i+1;
        plot(c.px(1:i),c.py(1:i),'k.');
        %plot(c{1}.before.px(i),c{1}.before.py(i),'k.');
        xlim([0 120]);ylim([0 120]);
        text(110,110,num2str(i),'fontsize',14,'color','k');
        text(110,110,num2str(i),'fontsize',14,'color','w');
        drawnow();
        pause(1/20000);
    end; toc
end















