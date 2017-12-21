function TimeTab(fig, tab, m)
    
    gui.Window = fig;
    gui.m = m;
    gui.m.tab = tab;
    gui.m.lag = 2.4;
    gui.m.binspike = 0.06;
    gui.m.sigma = 2;
    gui.m.version = 'n';
    gui.m.bindeg = 360;
    gui.m.sesh = 'before';
    gui.m.overlap = 0.002;
    gui.m.hamwin = 15;
    gui.Window.UserData.gidsFns{end+1} = @updateGids;
    %top box
    TimeTabBox = uix.HBoxFlex( 'Parent', tab, 'Spacing', 3);
    gui.paramPanel = uix.Panel('Parent', TimeTabBox,'Title', 'Parameters' );
    gui.viewPanel = uix.Panel('Parent', TimeTabBox,'Title', 'View','fontweight','bold',...
        'FontSize', 11, 'TitlePosition','centerTop','ForegroundColor','r');
    set( TimeTabBox, 'Widths', [200,-1] );% Adjust the main layout
    %% VIEW PANEL
    viewLayout = uix.VBox( 'Parent', gui.viewPanel,'Padding', 0, 'Spacing', 0);
    t = [];
    gui.viewContainer = uicontainer('Parent', viewLayout); t = [t -1]; %'backgroundColor','w'
    set(viewLayout, 'Heights', t);
    %% PARAMETERS PANEL
    paramsLayout = uix.VBox( 'Parent', gui.paramPanel,'Padding', 3, 'Spacing', 3 );
    %DEG
    degBoxV = uix.VBox( 'Parent', paramsLayout,'Padding', 3, 'Spacing', 3 );
    uicontrol('Parent', degBoxV,'Style','text',...
        'HorizontalAlignment', 'left','String', 'Direction Bins(deg):');
    degBoxBoxH = uix.HBox( 'Parent', degBoxV,'Padding', 3, 'Spacing', 3 );
    step = 3;step = step/(360-3); %max-min
    gui.degSlider = uicontrol( 'Style', 'slider', 'Min', 3,'Max', 360,'Parent', degBoxBoxH, ...
        'Value', 360,'SliderStep',[step 3*step],'Callback', @onDegSlide);
    gui.degEdit = uicontrol( 'Parent', degBoxBoxH, 'Style', 'edit', 'String', 360,...
        'Callback', @onDegEdit);
    addlistener(gui.degSlider,'ContinuousValueChange',@(hObject, event) onDegSliding(hObject, event));
    set( degBoxBoxH, 'Widths', [-4 -1] );
    %LAG
    lagBoxV = uix.VBox( 'Parent', paramsLayout,'Padding', 3, 'Spacing', 3 , 'Visible', 'on');
    uicontrol('Parent', lagBoxV,'Style','text',...
        'HorizontalAlignment', 'left','String', 'Correlation Lag(s):');
    lagBoxH = uix.HBox( 'Parent', lagBoxV,'Padding', 3, 'Spacing', 3 );
    step = 0.1; step = step/(5-0.1); %max-min
    gui.lagSlider = uicontrol( 'Style', 'slider', 'Min', 0.1,'Max', 5.0, ...
        'Parent', lagBoxH,'Value', 2.4,'SliderStep',[step 3*step],'Callback', @onLagSlide);
    gui.lagEdit = uicontrol( 'Parent', lagBoxH, 'Style', 'edit', 'String', 2.4,...
        'Callback', @onLagEdit);
    addlistener(gui.lagSlider,'ContinuousValueChange',@(hObject, event) onLagSliding(hObject, event));
    set( lagBoxH, 'Widths', [-4 -1] );
    %SPIKE
    spikeBinBoxV = uix.VBox( 'Parent', paramsLayout,'Padding', 3, 'Spacing', 3, 'Visible', 'on');
    uicontrol('Parent', spikeBinBoxV,'Style','text',...
        'HorizontalAlignment', 'left','String', 'Spike Time Bin(s):');
    spikeBinBoxH = uix.HBox( 'Parent', spikeBinBoxV, ...
        'Padding', 3, 'Spacing', 3 );
    step = 0.02; step = step/(0.2-0.02); %max-min
    gui.spikeBinSlider = uicontrol( 'Style', 'slider','Min', 0.02,'Max', 0.2, ...
        'Parent', spikeBinBoxH,'Value', 0.06,'SliderStep',[step 3*step],'Callback', @onSpikeBinSlide);
    gui.spikeBinEdit = uicontrol( 'Parent', spikeBinBoxH, 'Style', 'edit', 'String', 0.06, ...
        'Callback', @onSpikeBinEdit);
    addlistener(gui.spikeBinSlider,'ContinuousValueChange',@(hObject, event) onSpikeSliding(hObject, event));
    set( spikeBinBoxH, 'Widths', [-4 -1] );
    %SIGMA
    sigmaBoxV = uix.VBox( 'Parent', paramsLayout,'Padding', 3, 'Spacing', 3, 'Visible', 'on' );
    uicontrol('Parent', sigmaBoxV,'Style','text',...
        'HorizontalAlignment', 'left','String', 'Smoothing (Butterworth):');
    sigmaBoxH = uix.HBox( 'Parent', sigmaBoxV,'Padding', 3, 'Spacing', 3 );
    maxstep = 100;
    step = 1; step = step/(maxstep); %max-min
    gui.sigmaSlider = uicontrol( 'Style', 'slider','Min', 0,'Max', maxstep, 'SliderStep',[step 2*step],...
        'Parent', sigmaBoxH,'Value', gui.m.sigma,'Callback', @onSigmaSlide);
    gui.sigmaEdit = uicontrol( 'Parent', sigmaBoxH, 'Style', 'edit', 'String', gui.m.sigma, ...
        'Callback', @onSigmaEdit);
    addlistener(gui.sigmaSlider,'ContinuousValueChange',@(hObject, event) onSigmaSliding(hObject, event));
    set( sigmaBoxH, 'Widths', [-4 -1] );
    %OVERLAP
        overlapBoxV = uix.VBox( 'Parent', paramsLayout,'Padding', 3, 'Spacing', 3, 'Visible', 'on' );
    uicontrol('Parent', overlapBoxV,'Style','text',...
        'HorizontalAlignment', 'left','String', 'Remove overlapping spikes [0-0.5s]:');
    overlapBoxH = uix.HBox( 'Parent', overlapBoxV,'Padding', 3, 'Spacing', 3 );
    mnx = [0 0.5]; step = 0.001; step = step/(mnx(2)-mnx(1)); %max-min
    gui.overlapSlider = uicontrol( 'Style', 'slider','Min', mnx(1),'Max', mnx(2), 'SliderStep',[step 2*step],...
        'Parent', overlapBoxH,'Value', 0,'Callback', @onOverlapSlide);
    gui.overlapEdit = uicontrol( 'Parent', overlapBoxH, 'Style', 'edit', 'String', 0, ...
        'Callback', @onOverlapEdit);
    addlistener(gui.overlapSlider,'ContinuousValueChange',@(hObject, event) onOverlapSliding(hObject, event));
    set( overlapBoxH, 'Widths', [-4 -1] );   
    %HAMWIN
        hamwinBoxV = uix.VBox( 'Parent', paramsLayout,'Padding', 3, 'Spacing', 3, 'Visible', 'on' );
    uicontrol('Parent', hamwinBoxV,'Style','text',...
        'HorizontalAlignment', 'left','String', 'size of ham window:');
    hamwinBoxH = uix.HBox( 'Parent', hamwinBoxV,'Padding', 3, 'Spacing', 3 );
    mnx = [0 500]; step = 1; step = step/(mnx(2)-mnx(1)); %max-min
    gui.hamwinSlider = uicontrol( 'Style', 'slider','Min', mnx(1),'Max', mnx(2), 'SliderStep',[step 2*step],...
        'Parent', hamwinBoxH,'Value', 0,'Callback', @onHamwinSlide);
    gui.hamwinEdit = uicontrol( 'Parent', hamwinBoxH, 'Style', 'edit', 'String', 0, ...
        'Callback', @onHamwinEdit);
    addlistener(gui.hamwinSlider,'ContinuousValueChange',@(hObject, event) onHamwinSliding(hObject, event));
    set( hamwinBoxH, 'Widths', [-4 -1] );   
    %GROUP
    groupBoxV = uix.VBox( 'Parent', paramsLayout,'Padding', 3, 'Spacing', 3 );
    uicontrol('Parent', groupBoxV,'Style','text',...
        'HorizontalAlignment', 'left','String', 'Group:');
    gui.groupListBox = uicontrol( 'Style', 'list','BackgroundColor', 'w','Parent', groupBoxV, ...
        'String', gui.Window.UserData.gids,'Value', 1,'Callback', @onGroupListSelection);
    %MID
    midBoxV = uix.VBox( 'Parent', paramsLayout,'Padding', 3, 'Spacing', 3 );
    uicontrol('Parent', midBoxV,'Style','text','HorizontalAlignment', ...
        'left', 'String','Muscimol Bin:');
    gui.midListBox = uicontrol( 'Style', 'list','BackgroundColor', 'w', ...
        'Parent', midBoxV,'String', {'before';'midall';'after'}, ...
        'Value', 1,'Callback', @onMidListSelection, 'Enable', 'on');
    if ~gui.m.g(1).after.exists
        set(gui.midListBox,'String', {'before';'midall'});
    end
    %RUN
    runBoxV = uix.VBox( 'Parent', paramsLayout,'Padding', 3, 'Spacing', 3 );
    gui.RunButton = uicontrol( 'Style', 'PushButton','Parent', runBoxV,'String', 'Analysis!', ...
        'Callback', @onRunButton );
    gui.status = uicontrol('Parent', runBoxV,'Style','text', 'Tag', 'status',...
        'HorizontalAlignment', 'left', 'String', 'ok', 'FontSize',9);
    %% PARAMETERS PANEL TIGHTEN UP
    %Tighten UP
    t = []; a = 50;
    set(degBoxV, 'Heights', [20 -1]); t = [t a]; % Make the list fill the space
    set(lagBoxV, 'Heights', [20 -1]); t = [t a];% Make the list fill the space
    set(spikeBinBoxV, 'Heights', [20 -1]); t = [t a]; % Make the list fill the space
    set(sigmaBoxV, 'Heights', [20 -1]); t = [t a]; % Make the list fill the space
    set(overlapBoxV, 'Heights', [20 -1]); t = [t a]; % Make the list fill the space
    set(hamwinBoxV, 'Heights', [20 -1]); t = [t a]; % Make the list fill the space
    set(groupBoxV, 'Heights', [20 -1]); t = [t 200]; % Make the list fill the space
    set(midBoxV, 'Heights', [20 -1]); t = [t 120]; % Make the list fill the space
    set(runBoxV, 'Heights', [30 -1]); t = [t -1]; % Make the list fill the space
    
    set(paramsLayout, 'Heights', t); % Make the lists fill the space
    
    function updateGids()
        set(gui.groupListBox, 'String', gui.Window.UserData.gids);
    end
    
    %% CALLBACKS
    %DEG
    function onDegSlide( src, ~ )
        t = round2r(get(src, 'Value'),3);
        gui.m.bindeg = t;
        set(gui.degEdit, 'String', t);
        run();
    end
    function onDegEdit( src, ~ )
        t = round2r(str2double(get(src, 'String')),3);
        gui.m.bindeg = t;
        set(gui.degSlider, 'Value', t);
    end
    function onDegSliding(hObject,event)
        t = round2r(get(hObject,'Value'),3);
        set(gui.degEdit, 'String', t);
    end
    %LAG
    function onLagSlide( src, ~ )
        t = round2r(get(src, 'Value'), 0.1);
        gui.m.lag = t;
        set(gui.lagEdit, 'String', t);
        %run();
    end
    function onLagEdit( src, ~ )
        t = round2r(str2double(get(src, 'String')),0.1);
        gui.m.lag = t;
        set(gui.lagSlider, 'Value', t);
    end
    function onLagSliding(hObject,event)
        t = round2r(get(hObject,'Value'),0.1);
        set(gui.lagEdit, 'String', t);
    end
    %SPIKE
    function onSpikeBinSlide( src, ~ )
        t = round2r(get( src, 'Value'),0.02);
        gui.m.binspike = t;
        set(gui.spikeBinEdit, 'String', t);
        %run();
    end
    function onSpikeBinEdit(src, ~ )
        t = round2r(str2double(get(src, 'String')),0.02);
        gui.m.binspike =t;
        set(gui.spikeBinSlider, 'Value', t);
    end
    function onSpikeSliding(hObject,event)
        t = round(get(hObject,'Value')/0.02)*0.02;
        set(gui.spikeBinEdit, 'String', t);
    end
    %SIGMA
    function onSigmaSlide( src, ~ )
        t = round(get( src, 'Value'));
        gui.m.sigma = t;
        set(gui.sigmaEdit, 'String', t);
        %run();
    end
    function onSigmaEdit( src, ~ )
        t = round(str2double(get(src, 'String')));
        gui.m.sigma =t;
        set(gui.sigmaSlider, 'Value', t);
    end
    function onSigmaSliding(hObject,event)
        t = round2r(get(hObject,'Value'),1);
        set(gui.sigmaEdit, 'String', t);
    end
    %OVERLAP
    function onOverlapSlide( src, ~ )
        t = round(get( src, 'Value'),3);
        gui.m.overlap = t;
        set(gui.overlapEdit, 'String', t);
        %run();
    end
    function onOverlapEdit( src, ~ )
        t = round(str2double(get(src, 'String')),3);
        gui.m.overlap =t;
        set(gui.overlapSlider, 'Value', t);
    end
    function onOverlapSliding(hObject,event)
        t = round2r(get(hObject,'Value'),0.001);
        set(gui.overlapEdit, 'String', t);
    end
 %HAMWIN
    function onHamwinSlide( src, ~ )
        t = round(get( src, 'Value'));
        gui.m.hamwin = t;
        set(gui.hamwinEdit, 'String', t);
        %run();
    end
    function onHamwinEdit( src, ~ )
        t = round(str2double(get(src, 'String')),3);
        gui.m.hamwin =t;
        set(gui.hamwinSlider, 'Value', t);
    end
    function onHamwinSliding(hObject,event)
        t = round2r(get(hObject,'Value'),1);
        set(gui.hamwinEdit, 'String', t);
    end
    % GROUP FNS
    function onGroupListSelection( src, t )
        set(gui.status, 'ForegroundColor',[0,0,1]);
           t = cellstr(get( src, 'String' ));
        gui.m.gid = str2double(t(get( src, 'Value' )));
        gui.m.g = gui.m.groups{gui.m.gid};
        l = length(gui.m.g);
        set(gui.status, 'String', sprintf('group size: %d',l));
        set(gui.midListBox,'String',{'before';'midall';'after'},'Value', 1);
        %gui.m.mid = 1;
        %drawnow(); %TRY??
        run();
    end
    % MID
    function onMidListSelection( src, e )
        %gui.m.mid = get( src, 'Value' )-1;
        t = get( src, 'String');
        gui.m.sesh = str2double(t{get( src, 'Value' )});
        if isnan(gui.m.sesh )
            gui.m.sesh = t{get( src, 'Value' )};
        end
    end
    
    function onRunButton( ~, ~ )
        run();
    end % onDemoHelp
    function onExit( ~, ~ )
        % User wants to quit out of the application
        delete( gui.Window );
    end % onExit
    %% RUN
    function run()
        set(gui.status, 'ForegroundColor',[0,0,1]);
        disp(['time run g', num2str(gui.m.gid)]);
        delete(findobj(gui.m.tab,'type','axes'));
        drawnow();
        gui.m.grid_thresh = gui.Window.UserData.gridThresh;
        params = gui.m;
        set(gui.status, 'String', 'computing....');
        %enable(handles, false)
        params.sesh = params.sesh;
        params.parent = gui.viewContainer;
        params.fig = gui.Window;
        params.good = gui.m.g(gui.Window.UserData.cids{gui.m.gid});
        [~,v] = plotByDirectionMainTimeTab(params);
        set(gui.status, 'String', sprintf('%s',v));
        set(gui.status, 'ForegroundColor',[0,0,1]);
        set(gui.status, 'String', v);
        set(gui.viewPanel,'Title', v);
        
        %enable(handles, true);
    end
end


