function varargout = reviewRegistration2(varargin)
% REVIEWREGISTRATION2 MATLAB code for reviewRegistration2.fig
%      REVIEWREGISTRATION2, by itself, creates a new REVIEWREGISTRATION2 or raises the existing
%      singleton*.
%
%      H = REVIEWREGISTRATION2 returns the handle to a new REVIEWREGISTRATION2 or the handle to
%      the existing singleton*.
%
%      REVIEWREGISTRATION2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in REVIEWREGISTRATION2.M with the given input arguments.
%
%      REVIEWREGISTRATION2('Property','Value',...) creates a new REVIEWREGISTRATION2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before reviewRegistration2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to reviewRegistration2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help reviewRegistration2

% Last Modified by GUIDE v2.5 21-Mar-2017 15:35:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @reviewRegistration2_OpeningFcn, ...
                   'gui_OutputFcn',  @reviewRegistration2_OutputFcn, ...
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

function [sessionData, dateList] = loadSessionList(mouse,location)
% load a list of all the session dates (and if possible the time chunks)

% get the name struct for this mouse
nS = getNameStruct(mouse,'*',location);

localizationFileNameField = 'refLocBaselineFileName';

% get a list of files
locFNameTemp = nS.(localizationFileNameField);
locFNameTempStart = ...
    locFNameTemp(1:(strfind(locFNameTemp,'*')-1));
locFNameTempEnd = ...
    locFNameTemp((strfind(locFNameTemp,'*')+1):end);

flist = strsplit(ls(locFNameTemp));

% get a list of all the dates that have a localization file
dateList = strrep(strrep(flist, locFNameTempStart, ''),...
    locFNameTempEnd,'');
dateList = dateList(cellfun(@(x) ~isempty(x),dateList));
sessionData = cell(size(dateList));

% get the number of clusters associated with each session date
for dd = 1:length(dateList)
    nS = getNameStruct(mouse,dateList{dd},location);
    clustFile = load(nS.clusterFileName,'clusterSpec');
    nClusters = size(clustFile.clusterSpec,2) ;
%     sessionList{dd} = cell(1,nClusters);
    for cc = 1:nClusters
        sessionData{dd}{cc} = struct();
    end
end


function [handles] = loadDataset(handles)

nS = getNameStruct(handles.mouse,...
    handles.dateList{handles.currentSessionInd},handles.location);

% should create the buttons here

% get the automatically aligned block center positions
if ~isfield(handles.sessionData{handles.currentSessionInd...
        }{handles.currentChunkInd},'autoPts') ...
    || isempty(handles.sessionData{handles.currentSessionInd...
        }{handles.currentChunkInd}.autoPts,'var')
    
    autoPts = load(nS.refLocBaselineFileName,'xyzrcoClusterPeaks_auto');
%     handles.sessionData{handles.currentSessionInd ...
%         }{handles.currentChunkInd}.xyzrcoClusterPeaks_auto = ...
%         autoPts.xyzrcoClusterPeaks_auto;
    
    for cc = 1:length(handles.sessionData{handles.currentSessionInd})
        handles.sessionData{handles.currentSessionInd ...
        }{handles.currentChunkInd}.autoPts = num2cell(squeeze(autoPts.xyzrcoClusterPeaks_auto(1:3,...
        handles.currentChunkInd,:)),1);
    end
end


function [] = updatePlot(handles)
    
    autoPts = handles.sessionData{handles.currentSessionInd ...
        }{handles.currentChunkInd}.autoPts;
    autoPts = cat(2,autoPts{:});
    plot(handles.imRefAx,autoPts(1,:),autoPts(2,:),'.');
    %handles.overlayAx
        
    


function [] = createBlockButtonGrid(handles)
% 

gridPos = get(handles.blockGridAx,'position');

bg = uibuttongroup('Visible','off',...
                  'Position',gridPos);
handles.bg=bg;
              
% Create three radio buttons in the button group.
nX = 6; nY = 6; gapX = .1; gapY = .1; bounds = [0 0 1 1];
% identify subplot size
spacingX = bounds(3)/(nX+gapX/(1+gapX));
spacingY = bounds(4)/(nY+gapY/(1+gapY));
sizeX = spacingX/(1+gapX);
sizeY = spacingY/(1+gapY);


for xx = 1:nX
    for yy = 1:nY
    buttonNum = 6*(yy-1)+xx;
    handles.blockButton(yy,xx) = uicontrol(handles.bg,'Style',...
        'pushbutton',...
        'String',num2str(buttonNum),...
        'HandleVisibility','off',...
        'Units','normalized',...
        'Value',0,...
        'Position',[bounds(1)+spacingX*(gapX + xx-1) ...
        bounds(2)+bounds(4)-spacingY*yy sizeX sizeY],...
        'Callback',@bselection);
    end
end
              
axis(handles.blockGridAx,'off')
% Make the uibuttongroup visible after creating child objects. 
handles.bg.Visible = 'on';

% launch gui when click a block
function bselection(bg,event)
cpselect(randn(10),randn(100));



% --- Executes just before reviewRegistration2 is made visible.
function reviewRegistration2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to reviewRegistration2 (see VARARGIN)

% Record mouse name
handles.mouse = varargin{1};

% Record recording site location
handles.location = varargin{2};

% might want to optionally provide a date and a stack

% get list of sessions
[handles.sessionData, handles.dateList] = ...
    loadSessionList(handles.mouse,handles.location);

set(handles.sessionMenu,'String',handles.dateList)

handles.nChunksInSession = cellfun(@(x) length(x), handles.sessionData);
handles.sessionStartChunk = [1 1+cumsum(handles.nChunksInSession)];
handles.nChunks = sum(handles.nChunksInSession);
% use the session list to set up session slider
set(handles.sessionChunkSliderPage,'Max',handles.nChunks);
set(handles.sessionChunkSliderPage,'Value',1);
set(handles.sessionChunkSliderPage,'Min',1);
if handles.nChunks > 1
    set(handles.sessionChunkSliderPage,'SliderStep',[1 1]/(handles.nChunks-1));
end

handles.currentSessionInd  = 1;
handles.currentChunkInd = 1;

% load the first dataset
handles = loadDataset(handles);

% update the plot
updatePlot(handles);

% Choose default command line output for reviewRegistration2
handles.output = hObject;


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes reviewRegistration2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% function sessionChunkToInd(sess, chunk, nChunksInSession)
%     cumsum(nChunksInSession)
% 
% function indToSessionChunk


% --- Outputs from this function are returned to the command line.
function varargout = reviewRegistration2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function sliderZ_Callback(hObject, eventdata, handles)
% hObject    handle to sliderZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function sliderZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in imgRefToggle.
function imgRefToggle_Callback(hObject, eventdata, handles)
% hObject    handle to imgRefToggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of imgRefToggle


% --- Executes on button press in rigidToggle.
function rigidToggle_Callback(hObject, eventdata, handles)
% hObject    handle to rigidToggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rigidToggle


% --- Executes on button press in meanMultiZToggle.
function meanMultiZToggle_Callback(hObject, eventdata, handles)
% hObject    handle to meanMultiZToggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of meanMultiZToggle



function editZ_Callback(hObject, eventdata, handles)
% hObject    handle to editZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editZ as text
%        str2double(get(hObject,'String')) returns contents of editZ as a double


% --- Executes during object creation, after setting all properties.
function editZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ctrlPtButton.
function ctrlPtButton_Callback(hObject, eventdata, handles)
% hObject    handle to ctrlPtButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on slider movement.
function sessionChunkSliderPage_Callback(hObject, eventdata, handles)
% hObject    handle to sessionChunkSliderPage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
newPage = round(get(hObject,'Value'));

sessionNum = find(handles.sessionStartChunk <= newPage,1,'last');
handles.currentSessionInd = sessionNum;
handles.currentChunkInd = newPage - handles.sessionStartChunk(sessionNum) + 1;
set(handles.sessionChunkSliderPage, 'Value', newPage);

% update the session and chunk popup menu
set(handles.sessionMenu, 'Value', sessionNum);
set(handles.chunkMenu, 'String',  ...
    strsplit(num2str(1:handles.nChunksInSession(sessionNum))));
set(handles.chunkMenu, 'Value', handles.currentChunkInd);

% load data
handles = loadDataset(handles);

% replot
updatePlot(handles)

% handles = loadClusterPeaks(hObject,handles);
% plotBlockPositions(hObject,handles)



guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function sessionChunkSliderPage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sessionChunkSliderPage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in sessionMenu.
function sessionMenu_Callback(hObject, eventdata, handles)
% hObject    handle to sessionMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns sessionMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from sessionMenu

sessionNum = (get(hObject,'Value'));

handles.currentSessionInd = sessionNum;
handles.currentChunkInd = 1;

% update the session and chunk popup menu
set(handles.sessionMenu, 'Value', sessionNum);
set(handles.chunkMenu, 'String',  ...
    strsplit(num2str(1:handles.nChunksInSession(sessionNum))));
set(handles.chunkMenu, 'Value', handles.currentChunkInd);

% move the slider to the appropriate place
set(handles.sessionChunkSliderPage,'Value',handles.sessionStartChunk(sessionNum));

% load data
handles = loadDataset(handles);

% replot
updatePlot(handles)

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function sessionMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sessionMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in chunkMenu.
function chunkMenu_Callback(hObject, eventdata, handles)
% hObject    handle to chunkMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns chunkMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from chunkMenu
chunkNum = (get(hObject,'Value'));

handles.currentChunkInd = chunkNum;

% update chunk popup menu
set(handles.chunkMenu, 'Value', handles.currentChunkInd);

% move the slider to the appropriate place
set(handles.sessionChunkSliderPage,'Value',...
    handles.sessionStartChunk(handles.currentSessionInd)+handles.currentChunkInd-1);

% load data
handles = loadDataset(handles);

% replot
updatePlot(handles)

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function chunkMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to chunkMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in loadAllButton.
function loadAllButton_Callback(hObject, eventdata, handles)
% hObject    handle to loadAllButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
