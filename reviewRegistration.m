function varargout = reviewRegistration(varargin)
% REVIEWREGISTRATION MATLAB code for reviewRegistration.fig
%      REVIEWREGISTRATION, by itself, creates a new REVIEWREGISTRATION or raises the existing
%      singleton*.
%
%      H = REVIEWREGISTRATION returns the handle to a new REVIEWREGISTRATION or the handle to
%      the existing singleton*.
%
%      REVIEWREGISTRATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in REVIEWREGISTRATION.M with the given input arguments.
%
%      REVIEWREGISTRATION('Property','Value',...) creates a new REVIEWREGISTRATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before reviewRegistration_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to reviewRegistration_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help reviewRegistration

% Last Modified by GUIDE v2.5 28-Feb-2017 19:24:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @reviewRegistration_OpeningFcn, ...
    'gui_OutputFcn',  @reviewRegistration_OutputFcn, ...
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


% --- Executes just before reviewRegistration is made visible.
function reviewRegistration_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to reviewRegistration (see VARARGIN)

% Record mouse name
handles.mouse = varargin{1};

%
handles.location = 'L01';
handles.nS = getNameStruct(varargin{1},'%s', handles.location);

handles.summaryFileList = strsplit(ls(sprintf(...
    handles.nS.referenceLocalizationFileName,'*')));
handles.dateList = strrep(strrep(handles.summaryFileList,...
    [PATH_DATA handles.mouse '/'],''),...
    '/L01_referenceLocalization/summary.mat','');
stackPath = load(handles.summaryFileList{1},'stackPath');
handles.stackPath = fullfile(PATH_DATA,stackPath.stackPath);
stackIs   = imageSeries(handles.stackPath);
handles.stack = squeeze(permute(stackIs.images,[1 2 4 3]));
handles.roiSum = zeros(size(handles.stack(:,:,1,1)));
handles.bLine = zeros(size(handles.stack(:,:,1,1)));

handles.mouseDates = getMouseDates(handles.mouse);
handles.nPages     = length(handles.summaryFileList);
handles.nSlices    = size(handles.stack,3);

handles.currentDateInd = 1;
handles.currentDateStr = handles.dateList{handles.currentDateInd};
handles.currentZ = 25;

handles.currentCluster = 1;
handles.currentBlocks  = 1;

handles.currSumm = load(handles.summaryFileList{handles.currentDateInd});

handles.reregister=1;

handles.clusterFile = load(sprintf(handles.nS.clusterFileName,...
    handles.currentDateStr),'clusterSpec');
handles.nClusts = size(handles.clusterFile.clusterSpec,2);

handles.clustInfo = load(sprintf(handles.nS.clusterInfoFileNameFcn(...
    handles.currentCluster),handles.currentDateStr));
handles.nBlocks = size(handles.clustInfo.clusterBlockLocations,2);



% set up z slider
set(handles.sliderZ,'Max',handles.nSlices)
set(handles.sliderZ,'Value',handles.currentZ)
set(handles.sliderZ,'Min',1)
if handles.nSlices > 1
    set(handles.sliderZ,'SliderStep',[1 1]/(handles.nSlices-1))
end

% set up editable text field
set(handles.editZ,'String',num2str(handles.currentZ))

% set up date slider
set(handles.sliderPage,'Max',handles.nPages)
set(handles.sliderPage,'Value',handles.currentDateInd)
set(handles.sliderPage,'Min',1)
if handles.nPages > 1
    set(handles.sliderPage,'SliderStep',[1 1]/(handles.nPages-1))
end

% set up editable text field
set(handles.editPage,'String',num2str(handles.currentDateInd))

% Choose default command line output for reviewRegistration
handles.output = hObject;



% UIWAIT makes reviewRegistration wait for user response (see UIRESUME)
% uiwait(handles.figure1);
handles = loadClusterPeaks(hObject,handles);
plotBlockPositions(hObject,handles)
updatePlot(hObject, handles)

% Update handles structure
guidata(hObject, handles);

% give control to the slider
uicontrol(handles.sliderPage)




% --- Executes during object creation, after setting all properties.
function editPage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Outputs from this function are returned to the command line.
function varargout = reviewRegistration_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on slider movement.
function sliderPage_Callback(hObject, eventdata, handles)
% hObject    handle to sliderPage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
newPage = round(get(hObject,'Value'));

set(handles.editPage,'String',num2str(newPage))

handles.currentDateInd = newPage;
handles.currentDateStr = handles.dateList{handles.currentDateInd};


handles = loadClusterPeaks(hObject,handles);
plotBlockPositions(hObject,handles)

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function sliderPage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderPage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end





% --- Executes on button press in resetButton.
function resetButton_Callback(hObject, eventdata, handles)
% hObject    handle to resetButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in refitButton.
function refitButton_Callback(hObject, eventdata, handles)
% hObject    handle to refitButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end






function editPage_Callback(hObject, eventdata, handles)
% hObject    handle to editPage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editPage as text
%        str2double(get(hObject,'String')) returns contents of editPage as a double

% get requested frame
currentDateInd = str2double(get(hObject,'String'));

% if user input value is outside acceptable range
if isnan(currentDateInd) || currentDateInd < 1 || currentDateInd > handles.nPages
    % reset to previous value and don't update plot
    set(hObject,'String',num2str(get(handles.sliderPage,'Value')))
    return
end

% update slider
set(handles.sliderPage,'Value',currentDateInd)

% update handles struct
handles.currentDateInd = currentDateInd;
handles.currentDateStr = handles.dateList{handles.currentDateInd};
handles.clustInfo = load(sprintf(handles.nS.clusterInfoFileNameFcn(...
    handles.currentCluster),handles.currentDateStr));

updatePlot(hObject, handles);

% return control to the slider
uicontrol(handles.sliderPage)




% --- Executes on slider movement.
function sliderZ_Callback(hObject, eventdata, handles)
% hObject    handle to sliderZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
newZ = round(get(hObject,'Value'));

set(handles.editZ,'String',num2str(newZ))

handles.currentZ = newZ;

updatePlot(hObject,handles)

% --- Executes during object creation, after setting all properties.
function sliderZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function editZ_Callback(hObject, eventdata, handles)
% hObject    handle to editZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editZ as text
%        str2double(get(hObject,'String')) returns contents of editZ as a double

% get requested frame
currentZ = str2double(get(hObject,'String'));

% if user input value is outside acceptable range
if isnan(currentZ) || currentZ < 1 || currentZ > handles.nSlices
    % reset to previous value and don't update plot
    set(hObject,'String',num2str(get(handles.sliderZ,'Value')))
    return
end

% update slider
set(handles.sliderZ,'Value',currentZ)

% update handles struct
handles.currentZ = currentZ;

updatePlot(hObject, handles);

% return control to the slider
uicontrol(handles.sliderPage)


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





function handles = loadClusterPeaks(hObject, handles)


filenameToSave = [sprintf(handles.nS.prefix,handles.currentDateStr)...
    '_referenceLocalizationBaseline.mat'];

if exist(filenameToSave,'file')
    temp = load(filenameToSave,'xyzrcoClusterPeaks_auto');
    handles.xyzrcoClusterPeaks_auto = temp.xyzrcoClusterPeaks_auto;
    handles.xyzrcoClusterPeaks = temp.xyzrcoClusterPeaks_auto;
end





function plotBlockPositions(hObject,handles)
cla(handles.blAx)
hold(handles.blAx,'on')
set(handles.blAx,'color',[1 1 1].*0)

if handles.showRois
    for bb = 1:handles.nBlocks
        rectWidth = diff(handles.clustInfo.clusterBlockLocations{bb}.blockRangeX)+1;
        rectHeight = diff(handles.clustInfo.clusterBlockLocations{bb}.blockRangeY)+1;
    
        cellfile=load(sprintf(handles.nS.cellFileNameFcn(handles.currentCluster,bb),...
        handles.currentDateStr),'rois');
        roi_w = nanmax(cellfile.rois,[],3);
        
        block = false(handles.clustInfo.clusterParameters.imageSize);
        xStart = handles.xyzrcoClusterPeaks(1,handles.currentCluster,bb) - (rectWidth-1)/2;
        yStart = handles.xyzrcoClusterPeaks(2,handles.currentCluster,bb) - (rectHeight-1)/2;
        xEnd = handles.xyzrcoClusterPeaks(1,handles.currentCluster,bb) + (rectWidth-1)/2;
        yEnd = handles.xyzrcoClusterPeaks(2,handles.currentCluster,bb) + (rectHeight-1)/2;
        block(round(yStart:yEnd),round(xStart:xEnd)) = true;
        

        pos = rotateAndSelectBlock(block, block, ...
            handles.xyzrcoClusterPeaks(4,handles.currentCluster,bb));
        
    end
else

thisCmap = colormapRedBlue;
thisCmap = thisCmap(linspace(1,size(colormapRedBlue,1),handles.nSlices),:);

for bb = 1:handles.nBlocks
    rectWidth = diff(handles.clustInfo.clusterBlockLocations{bb}.blockRangeX)+1;
    rectHeight = diff(handles.clustInfo.clusterBlockLocations{bb}.blockRangeY)+1;
    
    drawShadedRect(rectWidth, rectHeight, ...
        handles.xyzrcoClusterPeaks(1,handles.currentCluster,bb) ,...
        handles.xyzrcoClusterPeaks(2,handles.currentCluster,bb), ...
        handles.xyzrcoClusterPeaks(4,handles.currentCluster,bb),...
        thisCmap(handles.xyzrcoClusterPeaks(3,handles.currentCluster,bb),:),...
        handles.xyzrcoClusterPeaks(end,handles.currentCluster,bb),...
        handles.blAx)
    
end

axis(handles.blAx, [-10 520 -10 520])
end




function updatePlot(hObject, handles)

handles.currentSliceIm = handles.stack(:,:,handles.currentZ);

if 0
    handles.reregister=0;
    [roisPlot, xyzrcoClusterPeaks, roisRigid,params,extras] = registerRois('J117', '2015-12-06', 'L01',...
        'whichClusters',1,'whichBlocks',[],'normalizeStack',1,'normalizeBlock',1,'loadedStack',handles.stack);
    
    handles.roiSum = zeros(size(handles.currentSliceIm));
    handles.bLine = zeros(size(handles.currentSliceIm));
    
    for bb=1:36
        
        rotBaseline = normalizeToZeroOne(double(roisPlot(handles.currentCluster,bb).baseline));
        
        ws =roisPlot(handles.currentCluster,bb).w;
        ws = bsxfun(@rdivide,ws,max(max(ws,[],1),[],2));
        wMax = max(ws,[],3);
        
        for ii=1:numel(roisPlot(handles.currentCluster,bb).x)
            
            xInd = roisPlot(handles.currentCluster,bb).x(ii);
            yInd = roisPlot(handles.currentCluster,bb).y(ii);
            
            if yInd <= size(handles.roiSum,1) && xInd <= size(handles.roiSum,2)
                handles.roiSum(yInd,xInd) = max(handles.roiSum(yInd,xInd),wMax(ii));
                handles.bLine(yInd,xInd) = rotBaseline(ii);
            end
            
        end
        
    end
% clear previous plot
cla(handles.blAx)
axis(handles.blAx, 'off')
colormap(handles.blAx,bone)
imagesc(handles.bLine,'parent',handles.blAx);
end


cla(handles.refAx)


imagesc(handles.currentSliceIm,...
    'parent',handles.refAx);
axis(handles.refAx, 'off')
colormap(handles.refAx,bone)




% set title
%set(handles.textTitle,'String',sprintf('This is page number %d',handles.currentDateInd))

% save any current values
guidata(hObject, handles);


% --- Executes on button press in toggleShowRoi.
function toggleShowRoi_Callback(hObject, eventdata, handles)
% hObject    handle to toggleShowRoi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of toggleShowRoi
handles.showRois=get(hObject,'Value');
plotBlockPositions(hObject,handles)
guidata(hObject, handles);
