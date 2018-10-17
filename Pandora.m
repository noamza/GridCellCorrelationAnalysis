function Pandora
    %{ 
    TO DO:
     - reconstruct grid cell after smoothing
     - see if this explains clustering in muscimol grid cells or
     - or explained by speed
     - miniscope, probe 
    %}
    gui = createInterface(createData()); %#ok<*NASGU>
end

    function m = createData()
        fn = sprintf('C:\\Noam\\Data\\muscimol\\noam\\mini%dmin.mat', 45);
        fn = sprintf('C:\\Noam\\Data\\muscimol\\noam\\cells_%dmin_d_patchtraj_rayleigh',45);
        %fn = sprintf('C:\\Noam\\Data\\muscimol\\noam\\cells_%dmin_c_midall_gridscore.mat',45); ORIGINAL results
        disp(['loading....  ' fn]); %ascii 48
        tic; cells = load(fn); cells = cells.cells; toc;
        groups = findSimultaneouslyRecordedCells(cells);
        for i = 1:length(groups)
        end
        m.groups = groups;
        m.gid = 1; %initialize to group
        m.g = m.groups{m.gid};
        m.sesh = 'before';
        %m.grid_thresh = 0.5;
    end % createData
    
    function gui = createInterface(m)
        
        gui = struct();
        gui.m = m;
        gui.Window = figure( ...
            'Position', [10,20, 1250, 650],'Name', 'Griddy',...
           'NumberTitle', 'off','MenuBar', 'none' );%, ...
          %  'Toolbar', 'none',  'HandleVisibility', 'off' );
        %GRID THRESH
        gui.Window.UserData.gridThresh = 0.5;
        gui.Window.UserData.gridThreshMid = 0.3;
        gui.Window.UserData.gidsFns = {};
        Tabs = uix.TabPanel( 'Parent',gui.Window, 'tabwidth', 100);
        t = {};
        gui.GroupTab = uix.Panel('Parent', Tabs, 'Tag', 'group tab'); t{end+1} = 'Group';
        gui.TimeTab  = uix.Panel('Parent', Tabs, 'Tag', 'time tab');  t{end+1} = 'Time Correlation';
        gui.DriftTab = uix.Panel('Parent', Tabs, 'Tag', 'drift tab'); t{end+1} = 'Drift';
        gui.AniTab   = uix.Panel('Parent', Tabs, 'Tag', 'ani tab');   t{end+1} = 'Replay';
        gui.CellTab  = uix.Panel('Parent', Tabs, 'Tag', 'cell tab');  t{end+1} = 'Cell';
        gui.MiscTab  = uix.Panel('Parent', Tabs, 'Tag', 'agg tab');   t{end+1} = 'Aggregate';

                 Tabs.Selection = 2; % <<<< START UP TAB                

        TabGroup(  gui.Window, gui.GroupTab ,m);
        TabTime(   gui.Window, gui.TimeTab  ,m);
        TabDrift(  gui.Window, gui.DriftTab ,m);
        TabMotion( gui.Window, gui.AniTab   ,m);
        TabCell(   gui.Window, gui.CellTab  ,m);
        TabMisc(   gui.Window, gui.MiscTab  ,m);
        Tabs.TabTitles = t;
    end