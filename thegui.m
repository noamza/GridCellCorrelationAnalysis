function varargout = thegui(varargin)
    
    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1; %'gui_Name',       mfilename, ...
    gui_State = struct('gui_Name',       mfilename, ...
        'gui_Singleton',  gui_Singleton, ...
        'gui_OpeningFcn', @gui_OpeningFcn, ...
        'gui_OutputFcn',  @gui_OutputFcn, ...
        'gui_LayoutFcn',  [] , ...
        'gui_Callback',   []);
    if nargin && ischar(varargin{1})
        gui_State.gui_Callback = str2func(varargin{1});
    end
    
    if nargout
        [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
    else
        gui_mainfcn(gui_State, varargin{:});
    end
    % End initialization code - DO NOT EDIT
    
    % --- Executes just before gui is made visible.
function gui_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to gui (see VARARGIN)
    
    % Choose default command line output for gui
    
    
    
    handles.output = hObject;
    set(handles.figure1, 'Name', 'Grid Time Correlation')
    % Update handles structure
    
    %handles.tabManager = TabManager( hObject );
    
    guidata(hObject, handles);
    
    
    
    initialize_gui(hObject, handles, false);
    
    % UIWAIT makes gui wait for user response (see UIRESUME)
    % uiwait(handles.figure1);
    
    
    % --- Outputs from this function are returned to the command line.
function varargout = gui_OutputFcn(hObject, eventdata, handles)
    % varargout  cell array for returning output args (see VARARGOUT);
    % hObject    handle to figure
    % handles    structure with handles and user data (see GUIDATA)
    
    % Get default command line output from handles structure
    varargout{1} = handles.output;
    
    
    % --------------------------------------------------------------------
function initialize_gui(fig_handle, handles, isreset)
    % If the model field is present and the reset flag is false, it means
    % we are we are just re-initializing a GUI by calling it from the cmd line
    % while it is up. So, bail out as we dont want to reset the data.
    if isfield(handles, 'model') && ~isreset
        %         return;
    end
    
    handles.model.bindeg = 15;
    handles.model.lag  = 2.4;
    handles.model.binspike  = 0.06;
    handles.model.sigma = 2;
    handles.model.version = 'n';
    handles.model.grid_thresh = 0.5;
    handles.model.fig = handles.figure1;
    handles.model.mid = 1;
    handles.model.gid = 1;
    
    
    set(handles.bindeg, 'String', handles.model.bindeg);
    set(handles.lag,  'String', handles.model.lag);
    set(handles.binspike, 'String', handles.model.binspike);
    set(handles.sigma, 'String', handles.model.sigma);
    set(handles.status, 'String', 'on');
    
    % Update handles structure
    guidata(handles.figure1, handles);
    enable(handles, true);
    %guidata(hObject, handles);
    
    % --- Executes during object creation, after setting all properties.
function bindeg_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
function bindeg_Callback(hObject, eventdata, handles)
    % Hints: get(hObject,'String') returns contents of bindeg as text
    %        str2double(get(hObject,'String')) returns contents of bindeg as a double
    bindeg = str2double(get(hObject, 'String'));
    if isnan(bindeg)
        set(hObject, 'String', 0);
        errordlg('Input must be a number','Error');
    end
    handles.model.bindeg = bindeg;
    guidata(hObject,handles);
    run(handles);
    
function lag_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
function lag_Callback(hObject, eventdata, handles)
    lag = str2double(get(hObject, 'String'));
    if isnan(lag)
        set(hObject, 'String', 0);
        errordlg('Input must be a number','Error');
    end
    handles.model.lag = lag;
    guidata(hObject,handles);
    run(handles);
    
function binspike_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
function binspike_Callback(hObject, eventdata, handles)
    binspike = str2double(get(hObject, 'String'));
    if isnan(binspike)
        set(hObject, 'String', 0);
        errordlg('Input must be a number','Error');
    end
    binspike = round2r(binspike,0.02);
    handles.model.binspike = binspike;
    guidata(hObject,handles);
    run(handles);
    
function sigma_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
function sigma_Callback(hObject, eventdata, handles)
    sigma = str2double(get(hObject, 'String'));
    if isnan(sigma)
        set(hObject, 'String', 0);
        errordlg('Input must be a number','Error');
    end
    handles.model.sigma = sigma;
    guidata(hObject,handles);
    run(handles);
    
function mid_Callback(hObject, eventdata, handles)
    mid = round(str2double(get(hObject, 'String')));
    if isnan(mid)
        set(hObject, 'String', 0);
        errordlg('Input must be a number','Error');
    end
    handles.model.mid = mid;
    guidata(hObject,handles);
    run(handles);
    
    % --- Executes during object creation, after setting all properties.
function mid_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to mid (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called
    
    % Hint: edit controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    
    % --- Executes on button press in runb.
function runb_Callback(hObject, eventdata, handles)
    % hObject    handle to runb (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    run(handles);
    
    %mass = handles.model.density * handles.model.volume;
    %set(handles.mass, 'String', mass);
    
function enable(handles,bool)
    t = struct2cell(handles);
    for i = 1:length(t)
        if isa(t{i},'matlab.ui.control.UIControl')
            %disp('disable')
            if bool
                t{i}.Enable = 'on';
            else
                t{i}.Enable = 'off';
            end
        end
    end
    drawnow();
    
    
    
function run(handles)
    delete(findobj(0,'type','axes'));
    params = handles.model;
    g =  params.groups{params.gid};
    params.grid_thresh = 0.5
    mid = params.mid;
    good = [g(1)]; goodmid = [g(1)];
    for j = 1:length(g)
        if g(j).before.gridscore > params.grid_thresh;
            good(end+1) = g(j);
        end;
        if true %g(j).middle{mid}.gridscore > 0.01
            %goodmid(end+1) = g(j);
        end
    end
    good = good(2:end); goodmid = goodmid(2:end);
    
    v = 'before';
    if ~isempty(good) && length(good) >= 2
        set(handles.status, 'String', 'computing....');
        enable(handles, false)
        params.sesh = v;
        params.parent = handles.uipanel2;
        params.parent.Title = v;
        [~,v] = plotByDirectionMain(params);
        set(handles.status, 'String', sprintf('%s',v));
        v = sprintf('mid%d',mid);
        params.sesh = mid;
        params.parent = handles.uipanel3;
        params.parent.Title = v;
        [~,v] = plotByDirectionMain(params);
        set(handles.status, 'ForegroundColor',[0,0,1]);
        set(handles.status, 'String', v);
    else
        set(handles.status, 'ForegroundColor',[1,0,0]);
        set(handles.status, 'String', sprintf('good cells[%s,%d]',v,length(good)));
    end
    enable(handles, true);
    
    
    % --- Executes on selection change in
function group_Callback(hObject, eventdata, handles)
    % Hints: contents = cellstr(get(hObject,'String')) returns groupDEAD contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from groupDEAD
    contents = cellstr(get(hObject,'String')); %returns group contents as cell array
    gid = str2num(contents{get(hObject,'Value')})% returns selected item from group
    handles.model.gid = gid;
    l = length(handles.model.groups{gid});
    set(handles.status, 'String', sprintf('group size: %d',l));
    guidata(hObject,handles)
    g = handles.model.groups{handles.model.gid};
    set(handles.midSL,'Max',length(g(1).middle));
    set(handles.midSL,'Value', 1); step = 1;
    step = step/(get(handles.midSL,'max')-get(handles.midSL,'Min'));
    set(handles.midSL,'SliderStep',[step 2*step]);
    set(handles.mid,'String',1);
    drawnow();
    run(handles);
    
    
function lame(hObject, eventdata, handles)
    disp('ugh GUIDE');
    
    % --- Executes during object creation, after setting all properties.
function group_CreateFcn(hObject, eventdata, handles)
    fn = sprintf('C:\\Noam\\Data\\muscimol\\noam\\cells_%dmin_b_midscorrect.mat', 10);
    fprintf('loading %s',fn); %ascii 48
    tic; cells = load(fn); cells = cells.cells; toc;
    groups = findSimultaneouslyRecordedCells(cells);
    handles.model.groups = groups;
    %k = fieldnames(handles.model.groups); %dt = 0.02;
    %handles.model.k = k;
    set(hObject,'String',cellstr(num2str((1:length(groups))')));
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    handles.model.gid = 1;
    guidata(hObject,handles);
    
    
    % --- Executes on slider movement.
function spikebinSL_Callback(hObject, eventdata, handles)
    % hObject    handle to spikebinSL (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    % Hints: get(hObject,'Value') returns position of slider
    %        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
    binspike = get(hObject,'Value');
    if isnan(binspike)
        set(hObject, 'String', 0);
        errordlg('Input must be a number','Error');
    end
    binspike = round2r(binspike, 0.02); %nearest to dt
    handles.model.binspike = binspike;
    set(handles.binspike, 'String', num2str(binspike));
    guidata(hObject,handles)
    run(handles)
    
    % --- Executes during object creation, after setting all properties.
function spikebinSL_CreateFcn(hObject, eventdata, handles)
    
    % Hint: slider controls usually have a light gray background.
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end
    set(hObject,'Min',0.02);set(hObject,'Max',0.2);
    set(hObject,'Value', 0.06);
    step = 0.02/(get(hObject,'max')-get(hObject,'Min'));
    set(hObject,'SliderStep',[step 2*step]);
    
    
    % --- Executes on slider movement.
function sigmaSL_Callback(hObject, eventdata, handles)
    sigma = get(hObject,'Value')
    sigma = round(sigma);
    handles.model.sigma = sigma;
    set(handles.sigma, 'String', num2str(sigma));
    guidata(hObject,handles)
    run(handles)
    
    % --- Executes during object creation, after setting all properties.
function sigmaSL_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to sigmaSL (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called
    
    % Hint: slider controls usually have a light gray background.
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end
    set(hObject,'Min',1);set(hObject,'Max',6);
    set(hObject,'Value', 2);
    step = 1/(get(hObject,'max')-get(hObject,'Min'));
    set(hObject,'SliderStep',[step 2*step]);
    
    
    % --- Executes on slider movement.
function lagSL_Callback(hObject, eventdata, handles)
    lag = get(hObject,'Value')
    lag = round2r(lag,0.1); %nearest to dt
    handles.model.lag = lag;
    set(handles.lag, 'String', num2str(lag));
    guidata(hObject,handles)
    run(handles)
    
    % --- Executes during object creation, after setting all properties.
function lagSL_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to lagSL (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called
    
    % Hint: slider controls usually have a light gray background.
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end
    set(hObject,'Min',0.1);set(hObject,'Max',5);
    set(hObject,'Value', 2.4);
    step = 0.1/(get(hObject,'max')-get(hObject,'Min'));
    set(hObject,'SliderStep',[step 2*step]);
    
    
    % --- Executes on slider movement.
function bindegSL_Callback(hObject, eventdata, handles)
    bindeg = get(hObject,'Value');
    bindeg = round(bindeg);
    handles.model.bindeg = bindeg;
    set(handles.bindeg, 'String', num2str(bindeg));
    guidata(hObject,handles)
    run(handles)
    
    
    % --- Executes during object creation, after setting all properties.
function bindegSL_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to bindegSL (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called
    
    % Hint: slider controls usually have a light gray background.
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end
    set(hObject,'Min',5);set(hObject,'Max',180);
    set(hObject,'Value', 15);
    step = 3/(get(hObject,'max')-get(hObject,'Min'));
    set(hObject,'SliderStep',[step 2*step]);
    
    % --- Executes on slider movement.
function midSL_Callback(hObject, eventdata, handles)
    mid = get(hObject,'Value');
    mid = round(mid);
    handles.model.mid = mid;
    set(handles.mid, 'String', num2str(mid));
    guidata(hObject,handles)
    %run(handles)
    params = handles.model;
    if  length(findobj(0,'Parent',handles.uipanel2)) > 1%isfield(params,'gid')
        delete(findobj(0,'Parent',handles.uipanel3));
        set(handles.status, 'String', 'computing....');
        enable(handles, false)
        v = sprintf('mid%d',mid);
        params.sesh = mid;
        params.parent = handles.uipanel3;
        params.parent.Title = v;
        [~,v] = plotByDirectionMain(params);
        set(handles.status, 'String', v);
    else
        run(handles)
    end
    enable(handles, true);
    
    
    
    % --- Executes during object creation, after setting all properties.
function midSL_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to bindegSL (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called
    
    % Hint: slider controls usually have a light gray background.
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end
    set(hObject,'Min',1);set(hObject,'Max',3);
    set(hObject,'Value', 1); step = 1;
    step = step/(get(hObject,'max')-get(hObject,'Min'));
    set(hObject,'SliderStep',[step 2*step]);
