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
    step = 1; step = step/(10); %max-min
    gui.sigmaSlider = uicontrol( 'Style', 'slider','Min', 0,'Max', 10, 'SliderStep',[step 2*step],...
        'Parent', sigmaBoxH,'Value', 1,'Callback', @onSigmaSlide);
    gui.sigmaEdit = uicontrol( 'Parent', sigmaBoxH, 'Style', 'edit', 'String', 1, ...
        'Callback', @onSigmaEdit);
    addlistener(gui.sigmaSlider,'ContinuousValueChange',@(hObject, event) onSigmaSliding(hObject, event));
    set( sigmaBoxH, 'Widths', [-4 -1] );
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
        'Parent', midBoxV,'String', ['before';cellstr((num2str(1:length(gui.m.g(1).middle),'%d'))');'midall';'after'], ...
        'Value', 1,'Callback', @onMidListSelection, 'Enable', 'on');
    %RUN
    runBoxV = uix.VBox( 'Parent', paramsLayout,'Padding', 3, 'Spacing', 3 );
    gui.RunButton = uicontrol( 'Style', 'PushButton','Parent', runBoxV,'String', 'Analysis!', ...
        'Callback', @onRunButton );
    gui.status = uicontrol('Parent', runBoxV,'Style','text', 'Tag', 'status',...
        'HorizontalAlignment', 'left', 'String', 'ok', 'FontSize',9);
    %% PARAMETERS PANEL TIGHTEN UP
    %Tighten UP
    t = [];
    set(degBoxV, 'Heights', [20 -1]); t = [t 50]; % Make the list fill the space
    set(lagBoxV, 'Heights', [20 -1]); t = [t 50];% Make the list fill the space
    set(spikeBinBoxV, 'Heights', [20 -1]); t = [t 50]; % Make the list fill the space
    set(sigmaBoxV, 'Heights', [20 -1]); t = [t 50]; % Make the list fill the space
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
    % GROUP FNS
    function onGroupListSelection( src, t )
        set(gui.status, 'ForegroundColor',[0,0,1]);
           t = cellstr(get( src, 'String' ));
        gui.m.gid = str2double(t(get( src, 'Value' )));
        gui.m.g = gui.m.groups{gui.m.gid};
        l = length(gui.m.g);
        set(gui.status, 'String', sprintf('group size: %d',l));
        set(gui.midListBox,'String',['before';cellstr((num2str(1:length(gui.m.g(1).middle),'%d'))');'midall';'after'],'Value', 1);
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
        disp(['yo what dude ', num2str(gui.m.gid)]);
        %delete(findobj(gui.Window,'type','axes'));
        delete(findobj(gui.m.tab,'type','axes'));
        drawnow();
        gui.m.grid_thresh = gui.Window.UserData.gridThresh;
        params = gui.m;
        set(gui.status, 'String', 'computing....');
        %enable(handles, false)
        params.sesh = params.sesh;
        params.parent = gui.viewContainer;
        %params.parent.Title = v;
        params.fig = gui.Window;
        [~,v] = plotByDirectionMain(params);
        set(gui.status, 'String', sprintf('%s',v));
        %v = sprintf('mid%d',mid);
        %params.sesh = mid;
        %params.parent = gui.mContainer;
        %params.parent.Title = v;
        %[~,v] = plotByDirectionMain(params);
        set(gui.status, 'ForegroundColor',[0,0,1]);
        set(gui.status, 'String', v);
        set(gui.viewPanel,'Title', v);
        
        %enable(handles, true);
    end
end










%{


    % + Create the panels
    controlPanel = uix.Panel('Parent', TimeTabBox,'Title', 'Parameters' );
    beforeMidBox = uix.HBoxFlex( 'Parent', TimeTabBox, 'Spacing', 3);
    gui.bPanel = uix.Panel('Parent', beforeMidBox,'Title', 'Before');
    gui.mPanel = uix.Panel('Parent', beforeMidBox,'Title', 'Mid');
    set( beforeMidBox, 'Widths', [-1,-1] ); %SET WIDTH
    %UICONTAINER TO group axes with bars/legends
    gui.bContainer = uicontainer('Parent', gui.bPanel );
    gui.mContainer = uicontainer('Parent', gui.mPanel );
    set( TimeTabBox, 'Widths', [200,-1] );% + Adjust the main layout
    %CONTROL MASTER BOX
    paramsLayout = uix.VBox( 'Parent', controlPanel,'Padding', 3, 'Spacing', 3 );
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
    lagBoxV = uix.VBox( 'Parent', paramsLayout,'Padding', 3, 'Spacing', 3 );
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
    spikeBinBoxV = uix.VBox( 'Parent', paramsLayout,'Padding', 3, 'Spacing', 3 );
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
    sigmaBoxV = uix.VBox( 'Parent', paramsLayout,'Padding', 3, 'Spacing', 3 );
    uicontrol('Parent', sigmaBoxV,'Style','text',...
        'HorizontalAlignment', 'left','String', 'Gaussian(sigma):');
    sigmaBoxH = uix.HBox( 'Parent', sigmaBoxV,'Padding', 3, 'Spacing', 3 );
    step = 1; step = step/(7); %max-min
    gui.sigmaSlider = uicontrol( 'Style', 'slider','Min', 0,'Max', 7, 'SliderStep',[step 3*step],...
        'Parent', sigmaBoxH,'Value', 1,'Callback', @onSigmaSlide);
    gui.sigmaEdit = uicontrol( 'Parent', sigmaBoxH, 'Style', 'edit', 'String', 1, ...
        'Callback', @onSigmaEdit);
    addlistener(gui.sigmaSlider,'ContinuousValueChange',@(hObject, event) onSigmaSliding(hObject, event));
    set( sigmaBoxH, 'Widths', [-4 -1] );
    %GROUP
    groupBoxV = uix.VBox( 'Parent', paramsLayout,'Padding', 3, 'Spacing', 3 );
    uicontrol('Parent', groupBoxV,'Style','text',...
        'HorizontalAlignment', 'left','String', 'Group:');
    gui.groupListBox = uicontrol( 'Style', 'list','BackgroundColor', 'w','Parent', groupBoxV, ...
        'String', num2str((1:length(gui.m.groups))'),'Value', 1,'Callback', @onGroupListSelection);
    %MID
    midBoxV = uix.VBox( 'Parent', paramsLayout,'Padding', 3, 'Spacing', 3 );
    uicontrol('Parent', midBoxV,'Style','text','HorizontalAlignment', ...
        'left', 'String','Muscimol Bin:');
    gui.midListBox = uicontrol( 'Style', 'list','BackgroundColor', 'w', ...
        'Parent', midBoxV,'String', num2str((1:length(gui.m.g(1).middle))'), ...
        'Value', 1,'Callback', @onMidListSelection);
    %RUN
    runBoxV = uix.VBox( 'Parent', paramsLayout,'Padding', 3, 'Spacing', 3 );
    gui.RunButton = uicontrol( 'Style', 'PushButton','Parent', runBoxV,'String', 'Analysis!', ...
        'Callback', @onRunButton );
    gui.m.status = uicontrol('Parent', runBoxV,'Style','text', 'Tag', 'status',...
        'HorizontalAlignment', 'left', 'String', 'ok', 'FontSize',9);
    
    
    %Tighten UP
    t = [];
    set(degBoxV, 'Heights', [20 -1]); t = [t 50]; % Make the list fill the space
    set(lagBoxV, 'Heights', [20 -1]); t = [t 50];% Make the list fill the space
    set(spikeBinBoxV, 'Heights', [20 -1]); t = [t 50]; % Make the list fill the space
    set(sigmaBoxV, 'Heights', [20 -1]); t = [t 50]; % Make the list fill the space
    set(groupBoxV, 'Heights', [20 -1]); t = [t 200]; % Make the list fill the space
    set(midBoxV, 'Heights', [20 -1]); t = [t 100]; % Make the list fill the space
    set(runBoxV, 'Heights', [30 -1]); t = [t -1]; % Make the list fill the space
    set(paramsLayout, 'Heights', t); % Make the list fill the space
    % + Create the view
    %gui.Baxes = axes( 'Parent', gui.bContainer );
    %gui.Maxes = axes( 'Parent', gui.bContainer );
 % createInterface

%DEG
function onDegSlide( src, ~ )
    t = round(get(src, 'Value'));
    gui.m.bindeg = t;
    set(gui.degEdit, 'String', t);
    run();
end
function onDegEdit( src, ~ )
    t = round(str2double(get(src, 'String')));
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
    run();
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
    run();
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
%IGMA
function onSigmaSlide( src, ~ )
    t = round(get( src, 'Value'));
    gui.m.sigma = t;
    set(gui.sigmaEdit, 'String', t);
    run();
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
%GROUP
function onGroupListSelection( src, t )
    set(gui.m.status, 'ForegroundColor',[0,0,1]);
    gui.m.gid = get( src, 'Value' );
    disp(['g: ', num2str(gui.m.gid)])
    gui.m.g = gui.m.groups{gui.m.gid};
    l = length(gui.m.g);
    set(gui.m.status, 'String', sprintf('group size: %d',l));
    set(gui.midListBox,'String',num2str((1:length(gui.m.g(1).middle))'),'Value', 1);
    gui.m.mid = 1;
    %drawnow(); %TRY??
    run();
end

function onMidListSelection( src, t )
    gui.m.mid = get( src, 'Value' );
    run();
end

%-------------------------------------------------------------------------%
function onRunButton( ~, ~ )
    run();
end % onDemoHelp

%-------------------------------------------------------------------------%
function onExit( ~, ~ )
    % User wants to quit out of the application
    delete( gui.Window );
end % onExit

%}
