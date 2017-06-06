 % - Define dummy data: 11 time series.
 t       = 0 : 0.1 : 10 ;
 data    = 2 * repmat( sin(t).', 1,11 ) + rand( length(t), 11 ) ;
 nSeries = size( data, 2 ) ;
 % - Build figure.
 figure() ;  clf ;
 set( gcf, 'Color', 'White', 'Unit', 'Normalized', ...
    'Position', [0.1,0.1,0.6,0.6] ) ;
 % - Compute #rows/cols, dimensions, and positions of lower-left corners.
 nCol = 4 ;  nRow = ceil( nSeries / nCol ) ;
 rowH = 0.58 / nRow ;  colW = 0.7 / nCol ;
 colX = 0.06 + linspace( 0, 0.96, nCol+1 ) ;  colX = colX(1:end-1) ;
 rowY = 0.1 + linspace( 0.9, 0, nRow+1 ) ;  rowY = rowY(2:end) ;
 % - Build subplots axes and plot data.
 for dId = 1 : nSeries
    rowId = ceil( dId / nCol ) ;
    colId = dId - (rowId - 1) * nCol ;
    axes( 'Position', [colX(colId), rowY(rowId), colW, rowH] ) ;
    plot( t, data(:,dId), 'b' ) ;
    grid on ;
    xlabel( '\theta(t) [rad]' ) ;  ylabel( 'Anomaly [m]' ) ;
    title( sprintf( 'Time series %d', dId )) ;    
 end
 % - Build title axes and title.
 axes( 'Position', [0, 0.95, 1, 0.05] ) ;
 set( gca, 'Color', 'None', 'XColor', 'White', 'YColor', 'White' ) ;
 text( 0.5, 0, 'My Nice Title', 'FontSize', 14', 'FontWeight', 'Bold', ...
      'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;









function titl = plotGroupPrivate2(parent, g, thresh, iter, bs, midmax, splitPages)
    prnt = false;
    rmt = 1.5; % Rate Map Threshold Hz
    debug = '';
    %get dimensions of plot
    %assums longest middle vector contains times of others
    %WHICH CELLS TO DISPLAY
    gridThresh =  thresh;
    good = [g(1)]; bad = [g(1)]; %find low gridscore cells
    for j = 1:length(g) %put low grid scores at end
        if g(j).before.gridscore > gridThresh
            good(end+1) = g(j);
        else
            bad(end+1) = g(j);
        end
    end
    g = good(2:end); %[good(2:end) bad(2:end)]; Only show good sells
    if ~isempty(g)
        g = sortByMiddle(g); %align times of cells
        %TODO: move to end if 0.3 gridness or less for before
        for j = 1:length(g(1).middle) %used to align middles
            t(j) = g(1).middle{j}.pt(1); %CHECK TIMES ALL ALIGN??
        end
        maxmid = min(length(t), midmax); %find index for after
        g = sortByEllipse(g); %display by ellipse size;
    end
    %PLOT
    %HEIGHT
    r = 6; %default plot length how many cells per page, NOT rows, 2 rows per cell with trajectory
    if (mod(length(g), r) ~= 0 && length(g) > 2*r) || length(g) == r+1
        r = r+1;
    end
    r = 6;
    %m = length(g(1).middle); m NOW midmax
    if ~splitPages
        r = length(g);
    end
    %WIDTH
    c = 2*(2+midmax) + 1; % width of figure: (before + max middle + after) * 2, for ratemap + autocorr + mcross
    tot = 0; %total cells printed
    page = 0; %height = 0;
    while tot < length(g)
        page = page + 1;
        rr = r;
        if r > (length(g) - tot)*2  %*2
            rr = (length(g) - tot)*2; %*2
        end
        titl = sprintf('i%d_cells%d_%dmin_pg%dof%d',iter,length(g),bs,page,ceil(2*length(g)/r));
        %colormap(ax, 'jet');
        %%% plot each cell %%%
        z = 0; %row for this page
        while z < rr
            z = z+1;
            tot = tot + 1; %total cells printed
            r1 = g(tot);
            mcross = {};
            %%%plot rm
            ax = subplot(r, c, c*(z-1) + 1, 'Parent', parent);
            gc = '';
            if r1.before.gridscore < gridThresh
                gc = '*';
            end
            imagesc(ax, r1.before.rm), title(ax, sprintf('%sc%d: t(-1) %.0fHz',gc, r1.ind, r1.before.max_r));
            mcross{1} = r1.before.rm; ind = 1 + midmax;
            %plot middle rms
            for i = 1: min(length(r1.middle), midmax)
                if  goodCell(r1.middle{i}, rmt)
                    %not(r1.middle{i}.max_r < rmt || r1.middle{i}.max_r == 50)
                    ind = find(t == r1.middle{i}.pt(1)); %align middles;
                    ax = subplot(r, c, c*(z-1) + 1 + ind, 'Parent', parent);
                    imagesc(ax, r1.middle{i}.rm);
                    title(ax, sprintf('%.0f-%.0f(%0.fHz)', r1.middle{i}.pt(1)/60,...
                        r1.middle{i}.pt(length(r1.middle{i}.pt))/60, r1.middle{i}.max_r));
                    %title(ax, sprintf('%.0fHz', r1.middle{i}.max_r));
                    mcross{end + 1} = r1.middle{i}.rm;
                else
                    ind = find(t == r1.middle{i}.pt(1)); %align middles;
                    ax = subplot(r, c, c*(z-1) + 1 + ind,'Parent', parent);
                    %imagesc(ax, -1);
                    %title(ax, sprintf('%.0f-%.0f(>%0.fHz)', r1.middle{i}.pt(1)/60,...
                    %r1.middle{i}.pt(length(r1.middle{i}.pt))/60, rmt));
                end
                
            end                           %after
            if r1.after.exists == 1 %will be fixed
                %ax = subplot('Parent', parent, r, c, c*(z-1) + 1 + m + 1); %always plot at fixed
                ax = subplot(r, c, c*(z-1) + 1 + maxmid + 1,'Parent', parent); %plot after last mid
                imagesc(ax, r1.after.rm);title(ax, sprintf('t(+1) %.0fHz', r1.after.max_r));
                mcross{end + 1} = r1.after.rm;
            end
            %%%plot ac 
            plotAcorrModule(r1.before, parent, [r, c, c*(z-1) + 3 + midmax], 't(-1):');
            %plot middle ac
            for i = 1: min(length(r1.middle), midmax)
                if goodCell(r1.middle{i}, rmt)
                    ind = find(t == r1.middle{i}.pt(1));
                    %ax = subplot('Parent', parent, r, c, c*(z-1) + 1 + ind);
                    plotAcorrModule(r1.middle{i}, parent, [r, c, c*(z-1) + 3 + midmax + ind], '');
                else
                    ind = find(t == r1.middle{i}.pt(1)); %align middles;
                    ax = subplot(r, c, c*(z-1) + 3 + midmax + ind,'Parent', parent); %(m,n,l);
                    %title(ax, sprintf('%s%.2f(%d)','', r1.middle{i}.gridscore, length((r1.middle{i}.sx))));
                end
            end                           %after
            if r1.after.exists == 1  %will be fixed
                %ax = subplot('Parent', parent, r, c, c*(z-1) + 3 + 2*m + 1);
                %plotAcorrModule(r1.after, h, [r, c, c*(z-1) + 2 + 2*m + 2], 't(+1):'); %always plot last
                plotAcorrModule(r1.after, parent, [r, c, c*(z-1) + 3 + midmax + maxmid + 1], 't(+1):'); %plot after mid
                
            end
            %mcross
            ax = subplot(r, c, c*(z-1) + 2 + 2*midmax + 3,'Parent', parent); mc = allCorr(mcross,mcross);
            mc =  padarray(mc,midmax+2-size(mc),-1,'post'); 
            imagesc(ax, mc);  
            caxis(ax, [-1 1]);
            title(ax, sprintf('[%.1f %.1f]', min(mc(:)), max(mc(:))));
            
            %%%%%%%%%%%%%%% TRAJECTORY
            z = z+1;
            ax = subplot(r, c, c*(z-1) + 1,'Parent', parent);
            plot(ax, r1.before.px, flip(r1.before.py));
            hold(ax, 'on');
            scatter(ax, r1.before.sx, flip(r1.before.sy), '.'),...
                xlim(ax, [0 100]), ylim(ax, [0 100]);
            title(ax, sprintf('%.2f', r1.before.gridscore));
            mcrossT{1} = r1.before.rm;
            %plot middle rms
            for i = 1: min(length(r1.middle), midmax)
                if  goodCell(r1.middle{i}, rmt)
                    %not(r1.middle{i}.max_r < rmt || r1.middle{i}.max_r == 50)
                    ind = find(t == r1.middle{i}.pt(1)); %align middles;
                    ax = subplot(r, c, c*(z-1) + 1 + ind, 'Parent', parent);
                    plot(ax, r1.middle{i}.px, flip(r1.middle{i}.py));
                    hold(ax, 'on');
                    scatter(ax, r1.middle{i}.sx, flip(r1.middle{i}.sy), '.');
                    xlim(ax, [0 100]); ylim(ax, [0 100]);
                    title(ax, sprintf('%.2f', r1.middle{i}.gridscore));
                    %title(ax, sprintf('%.0fHz', r1.middle{i}.max_r));
                    mcrossT{end + 1} = r1.middle{i}.rm;
                else
                    ind = find(t == r1.middle{i}.pt(1)); %align middles;
                    ax = subplot(r, c, c*(z-1) + 1 + ind,'Parent', parent);
                end
                
            end                           %after
            if r1.after.exists == 1 %will be fixed
                %ax = subplot('Parent', parent, r, c, c*(z-1) + 1 + m + 1);
                ax = subplot(r, c, c*(z-1) + 1 + maxmid + 1,'Parent', parent); %plot after last mid
                plot(ax, r1.after.px, flip(r1.after.py));
                hold(ax, 'on');
                scatter(ax, r1.after.sx, flip(r1.after.sy), '.'),...
                    xlim(ax, [0 100]), ylim(ax, [0 100]);
                title(ax, sprintf('%.2f', r1.after.gridscore));
                mcrossT{end + 1} = r1.after.rm;
            end
            %time correlation
            %           ax = subplot('Parent', parent, r, c, c*(z-1) + 3 + m);
            %           t2 = g(j).before.st;
            %           [p, c] = time_correlation(t1,t2); 
            %MAKE AXIS NICE
            
            ax = findobj(parent,'Type','Axes');
            for i=1:length(ax)
                set(ax(i),'LooseInset', get(ax(i),'TightInset'));
                set(ax(i),'FontSize',7);
                axis(ax(i),'equal')
                axis(ax(i),'off')
                colormap(ax(i),'jet')
                %axis off; axis equal;
            end
        end %% END ROW %%
        %tightfig(parent);
        
        % END PLOTTING PAGE %%
        %[maxc maxl] = plotGroupTemporralCorrs(g, 500);
        %ax = subplot('Parent', parent, r,c,r*c);
        % *********
        sdir = 'C:\Noam\Output\muscimol\groups\';
        filename = sprintf('%s%s%s.png',sdir, debug, titl);
        %titl = sprintf('cells%d_pg%dof%d\n', length(g), page, ceil(length(g)/r));
        disp(filename);
        if prnt
            parent.PaperPositionMode = 'auto'; %fig %fig = gcf;
            print(filename, '-dpng','-r0');
            close(parent);
        end
        %}
    end %END PLOTTING GROUP
    %}
end



function demoBrowser()
%demoBrowser: an example of using layouts to build a user interface
%
%   demoBrowser() opens a simple GUI that allows several of MATLAB's
%   built-in demos to be viewed. It aims to demonstrate how multiple
%   layouts can be used to create a good-looking user interface that
%   retains the correct proportions when resized. It also shows how to
%   hook-up callbacks to interpret user interaction.
%
%   See also: <a href="matlab:doc Layouts">Layouts</a>

%   Copyright 2010-2013 The MathWorks, Inc.

% Data is shared between all child functions by declaring the variables
% here (they become global to the function). We keep things tidy by putting
% all GUI stuff in one structure and all data stuff in another. As the app
% grows, we might consider making these objects rather than structures.
data = createData();
gui = createInterface( data.DemoNames );

% Now update the GUI with the current data
updateInterface();
redrawDemo();

% Explicitly call the demo display so that it gets included if we deploy
displayEndOfDemoMessage('')

%-------------------------------------------------------------------------%
    function data = createData()
        % Create the shared data-structure for this application
        demoList = {
            'Complex surface'            'cplxdemo'
            'Cruller'                    'cruller'
            'Earth'                      'earthmap'
            'Four linked tori'           'tori4'
            'Klein bottle'               'xpklein'
            'Klein bottle (1)'           'klein1'
            'Knot'                       'knot'
            'Logo'                       'logo'
            'Spherical Surface Harmonic' 'spharm2'
            'Werner Boy''s Surface'      'wernerboy'
            };
        selectedDemo = 8;
        data = struct( ...
            'DemoNames', {demoList(:,1)'}, ...
            'DemoFunctions', {demoList(:,2)'}, ...
            'SelectedDemo', selectedDemo );
    end % createData

%-------------------------------------------------------------------------%
    function gui = createInterface( demoList )
        % Create the user interface for the application and return a
        % structure of handles for global use.
        gui = struct();
        % Open a window and add some menus
        gui.Window = figure( ...
            'Name', 'Gallery browser', ...
            'NumberTitle', 'off', ...
            'MenuBar', 'none', ...
            'Toolbar', 'none', ...
            'HandleVisibility', 'off' );
        
        % + File menu
        gui.FileMenu = uimenu( gui.Window, 'Label', 'File' );
        uimenu( gui.FileMenu, 'Label', 'Exit', 'Callback', @onExit );
        
        % + View menu
        gui.ViewMenu = uimenu( gui.Window, 'Label', 'View' );
        for ii=1:numel( demoList )
            uimenu( gui.ViewMenu, 'Label', demoList{ii}, 'Callback', @onMenuSelection );
        end
        
        % + Help menu
        helpMenu = uimenu( gui.Window, 'Label', 'Help' );
        uimenu( helpMenu, 'Label', 'Documentation', 'Callback', @onHelp );
        
        
        % Arrange the main interface
        mainLayout = uix.HBoxFlex( 'Parent', gui.Window, 'Spacing', 3 );
        
        % + Create the panels
        controlPanel = uix.BoxPanel( ...
            'Parent', mainLayout, ...
            'Title', 'Select a demo:' );
        gui.ViewPanel = uix.BoxPanel( ...
            'Parent', mainLayout, ...
            'Title', 'Viewing: ???', ...
            'HelpFcn', @onDemoHelp );
        gui.ViewContainer = uicontainer( ...
            'Parent', gui.ViewPanel );        

        % + Adjust the main layout
        set( mainLayout, 'Widths', [-1,-2]  );
               
        % + Create the controls
        controlLayout = uix.VBox( 'Parent', controlPanel, ...
            'Padding', 3, 'Spacing', 3 );
        gui.ListBox = uicontrol( 'Style', 'list', ...
            'BackgroundColor', 'w', ...
            'Parent', controlLayout, ...
            'String', demoList(:), ...
            'Value', 1, ...
            'Callback', @onListSelection);
        gui.HelpButton = uicontrol( 'Style', 'PushButton', ...
            'Parent', controlLayout, ...
            'String', 'Help for <demo>', ...
            'Callback', @onDemoHelp );
        set( controlLayout, 'Heights', [-1 28] ); % Make the list fill the space
        
        % + Create the view
        p = gui.ViewContainer;
        gui.ViewAxes = axes( 'Parent', p );
        
        
    end % createInterface

%-------------------------------------------------------------------------%
    function updateInterface()
        % Update various parts of the interface in response to the demo
        % being changed.
        
        % Update the list and menu to show the current demo
        set( gui.ListBox, 'Value', data.SelectedDemo );
        % Update the help button label
        demoName = data.DemoNames{ data.SelectedDemo };
        set( gui.HelpButton, 'String', ['Help for ',demoName] );
        % Update the view panel title
        set( gui.ViewPanel, 'Title', sprintf( 'Viewing: %s', demoName ) );
        % Untick all menus
        menus = get( gui.ViewMenu, 'Children' );
        set( menus, 'Checked', 'off' );
        % Use the name to work out which menu item should be ticked
        whichMenu = strcmpi( demoName, get( menus, 'Label' ) );
        set( menus(whichMenu), 'Checked', 'on' );
    end % updateInterface

%-------------------------------------------------------------------------%
    function redrawDemo()
        % Draw a demo into the axes provided
        
        % We first clear the existing axes ready to build a new one
        if ishandle( gui.ViewAxes )
            delete( gui.ViewAxes );
        end
        
        % Some demos create their own figure. Others don't.
        fcnName = data.DemoFunctions{data.SelectedDemo};
        switch upper( fcnName )
            case 'LOGO'
                % These demos open their own windows
                evalin( 'base', fcnName );
                gui.ViewAxes = gca();
                fig = gcf();
                set( fig, 'Visible', 'off' );
                
            otherwise
                % These demos need a window opening
                fig = figure( 'Visible', 'off' );
                evalin( 'base', fcnName );
                gui.ViewAxes = gca();
        end
        % Now copy the axes from the demo into our window and restore its
        % state.
        cmap = colormap( gui.ViewAxes );
        set( gui.ViewAxes, 'Parent', gui.ViewContainer );
        colormap( gui.ViewAxes, cmap );
        rotate3d( gui.ViewAxes, 'on' );
        % Get rid of the demo figure
        close( fig );
    end % redrawDemo

%-------------------------------------------------------------------------%
    function onListSelection( src, ~ )
        % User selected a demo from the list - update "data" and refresh
        data.SelectedDemo = get( src, 'Value' );
        updateInterface();
        redrawDemo();
    end % onListSelection

%-------------------------------------------------------------------------%
    function onMenuSelection( src, ~ )
        % User selected a demo from the menu - work out which one
        demoName = get( src, 'Label' );
        data.SelectedDemo = find( strcmpi( demoName, data.DemoNames ), 1, 'first' );
        updateInterface();
        redrawDemo();
    end % onMenuSelection


%-------------------------------------------------------------------------%
    function onHelp( ~, ~ )
        % User has asked for the documentation
        doc layout
    end % onHelp

%-------------------------------------------------------------------------%
    function onDemoHelp( ~, ~ )
        % User wnats documentation for the current demo
        showdemo( data.DemoFunctions{data.SelectedDemo} );
    end % onDemoHelp

%-------------------------------------------------------------------------%
    function onExit( ~, ~ )
        % User wants to quit out of the application
        delete( gui.Window );
    end % onExit

end % EOF


function start
set(handles.pushbutton4,'Visible','off');

% Choose default command line output for Projectile_UnitConv4
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Projectile_UnitConv4 wait for user response (see UIRESUME)
% uiwait(handles.figure1);

varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)

set(handles.pushbutton4,'Visible','on');

for k=1:length(tplot)
    plot(x(k),y(k),'*r','MarkerFaceColor','red','MarkerEdgeColor','k')
    plot(x(1),y(1),'*r','MarkerFaceColor','red');
    legend('Path of Object','Position of Object','Maximum Height',12);
    
    if y(k)==max(y)|| x(k)>=(max(x)/2)
        plot((max(x)/2),max(y),'*r');
    end

    pause(ptime);
    hold off
    
    drawnow
    if isappdata(handles.figure1,'stopPlot')
        break
    end
end
rmappdata(handles.figure1,'stopPlot')


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(handles.figure1,'stopPlot',1);
% Update handles structure
guidata(hObject, handles);