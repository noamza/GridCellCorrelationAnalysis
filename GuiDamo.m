function GuiDamo()

%GUI = Figure + data structure

gui = struct();

gui.window = figure(...
    'Position', [50,50, 500, 500],'Name', 'gui time !', ...
    'NumberTitle', 'off','MenuBar', 'none','Toolbar', 'none');

%global data
gui.win.UserData.global = 42;

tabs = uix.TabPanel( 'Parent',gui.window);

gui.plotTab  = uix.Panel('Parent', tabs);
gui.aniTab   = uix.Panel('Parent', tabs);

tabs.TabTitles = {'plot','replay'};

plotTab(gui.window, gui.plotTab);
aniTab(gui.window, gui.aniTab);

end

function plotTab(fig, tab)

tabdata = struct();
tabdata.t = 0:0.01:(2*pi);
tabdata.x = 1;

tabHBox = uix.HBoxFlex( 'Parent', tab, 'Padding', 2, 'Spacing',3);
controlPanel = uix.Panel('Parent', tabHBox,'Title', 'parameters');
%,'FontSize', 11, 'TitlePosition','centerTop');
viewPanel = uix.Panel('Parent', tabHBox,'Title', 'the view panel');
set( tabHBox, 'Widths', [200 -1] );

%% CONTROL PANEL
controlBoxV = uix.VBox( 'Parent', controlPanel, 'Spacing', 3);
%slider
slider = uicontrol(  'Parent',controlBoxV,'Style', 'slider', 'Min', 1,'Max', 10, 'Value',tabdata.x,...
    'Callback', @onSlide); %'SliderStep',[1/9 3/9],
%run button
runButton = uicontrol('Parent', controlBoxV,'Style', 'PushButton','String', 'Yala','Callback',@onRunButton);
%set heights of ui elements
set(controlBoxV,'Heights',[20,-1]);


%% CALLBACKS
    function onSlide(src,event)
        tabdata.x = get(src, 'Value')
    end

    function onRunButton( ~, ~ )
        run();
    end % onDemoHelp

    function run()
        delete(findobj(tab,'type','axes'));
        ax = axes('Parent',viewPanel);
        plot(ax,tabdata.t, sin(tabdata.x*tabdata.t));
        title(ax,sprintf('x = %.1f',tabdata.x));
        axis(ax,'square');axis(ax,'tight');
        
    end

end

function aniTab(fig, tab);
tabdata = struct();
tabdata.t = 0:0.01:(2*pi);
tabdata.x = 1;
tabdata.stop = false;

tabHBox = uix.HBoxFlex( 'Parent', tab, 'Padding', 2, 'Spacing',3);
controlPanel = uix.Panel('Parent', tabHBox,'Title', 'parameters');
%,'FontSize', 11, 'TitlePosition','centerTop');
viewPanel = uix.Panel('Parent', tabHBox,'Title', 'the view panel');
set( tabHBox, 'Widths', [200 -1] );

%% CONTROL PANEL
controlBoxV = uix.VBox( 'Parent', controlPanel, 'Spacing', 3);
%slider
slider = uicontrol(  'Parent',controlBoxV,'Style', 'slider', 'Min', 1,'Max', 10, 'Value',tabdata.x,...
    'Callback', @onSlide); %'SliderStep',[1/9 3/9],
%run button
uicontrol('Parent', controlBoxV,'Style', 'PushButton','String', 'Yala','Callback',@onRunButton);
uicontrol('Parent', controlBoxV,'Style', 'PushButton','String', 'Stop!','Callback',@onStopButton);
%set heights of ui elements
set(controlBoxV,'Heights',[20,-1,-1]);
    function onSlide(src,event)
        tabdata.x = get(src, 'Value')
    end

    function onRunButton( ~, ~ )
        run();
    end % onDemoHelp
    function onStopButton( ~, ~ )
        tabdata.stop = true;
    end % onDemoHelp
    
    function run(~,~)
        delete(findobj(tab,'type','axes'));
        ax = axes('Parent',viewPanel);
        tabdata.stop = false;
        %title(ax,sprintf('x = %.1f',tabdata.x));
        i = 1;
        while true
            i = i + 1
            if(i>length(tabdata.t))
                i = 1;
            end
            t = [tabdata.t(i+1:length(tabdata.t)), tabdata.t(1:i)];
            plot(ax,tabdata.t, sin(tabdata.x.*t));
            axis(ax,'square');axis(ax,'tight');
            drawnow;
            if tabdata.stop
                tabdata.stop = false;
                cla(ax);
                break
            end
        end
    end

end
