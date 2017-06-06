function Pandora
    
    gui = createInterface(createData());
    function m = createData()
        fn = sprintf('C:\\Noam\\Data\\muscimol\\noam\\mini%dmin.mat', 10);
        %fn = sprintf('C:\\Noam\\Data\\muscimol\\noam\\cells_%dmin_b_midscorrect.mat', 10);
        %fn = sprintf('C:\\Noam\\Data\\muscimol\\noam\\cells_%dmin_a.mat', 10);
        disp(['loading....  ' fn]); %ascii 48
        tic; cells = load(fn); cells = cells.cells; toc;
        groups = findSimultaneouslyRecordedCells(cells);
        m.groups = groups;
        m.gid = 1;
        m.g = m.groups{m.gid};
        m.sesh = 'before';
        %m.grid_thresh = 0.5;
    end % createData
    
    function gui = createInterface(m)        
        gui = struct();
        gui.m = m;
        gui.Window = figure( ...
            'Position', [100,100, 1050, 800],'Name', 'Griddy', ...
            'NumberTitle', 'off','MenuBar', 'none', ...
            'Toolbar', 'none',  'HandleVisibility', 'off' );
        gui.Window.UserData.gridThresh = 0.5;
        gui.Window.UserData.gidsFns = {};
        Tabs = uix.TabPanel( 'Parent',gui.Window, 'tabwidth', 100);
        t = {};
        gui.GroupTab   = uix.Panel('Parent', Tabs, 'Tag', 'group tab'); t{end+1} = 'Group';
        gui.TimeTab  = uix.Panel('Parent', Tabs, 'Tag', 'time tab'); t{end+1} = 'Time Correlation';
        gui.DriftTab = uix.Panel('Parent', Tabs, 'Tag', 'drift tab'); t{end+1} = 'Drift';
        gui.AniTab   = uix.Panel('Parent', Tabs, 'Tag', 'ani tab'); t{end+1} = 'Replay';

                 Tabs.Selection = 4; % <<<< START UP TAB                
        
        GroupTab(gui.Window, gui.GroupTab, m);
        TimeTab(gui.Window, gui.TimeTab, m);
        DriftTab(gui.Window, gui.DriftTab, m);
        MotionTab(gui.Window, gui.AniTab, m);
        Tabs.TabTitles = t;
        
        
        
    end
end