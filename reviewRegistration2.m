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

% Last Modified by GUIDE v2.5 11-Apr-2017 18:31:49

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


set(handles.loadingText, 'Visible', 'on');
drawnow
handles.nS = getNameStruct(handles.mouse,...
    handles.dateList{handles.currentSessionInd},handles.location);

sessInd = handles.currentSessionInd;
chunkInd = handles.currentChunkInd;


% load the stack
if ~isfield(handles, 'reference') || isempty(handles.reference)
    refDate = defaultStackDate(handles.mouse);
    refFilename = fullfile(PATH_DATA, handles.mouse, sprintf('reference_stack_%s.tif', refDate));
    refIs = imageSeries(refFilename);
    handles.reference = refIs.images;
end


if ~isfield(handles, 'chunkParams') || sessInd > length(handles.chunkParams) || ...
        chunkInd > length(handles.chunkParams{sessInd}) || ...
        isempty(handles.chunkParams{sessInd}{chunkInd})
    handles.chunkParams{sessInd}{chunkInd} = load(handles.nS.clusterInfoFileNameFcn(chunkInd));
end


%handles.reference = load(handles.chunkParams)

% should create the buttons here

% get the automatically aligned block center positions
if ~isfield(handles.sessionData{sessInd}{chunkInd},'autoPts') ...
    || isempty(handles.sessionData{sessInd}{chunkInd}.autoPts)
    
    autoPts = load(handles.nS.refLocBaselineFileName,'xyzrcoClusterPeaks_auto');
    
    for cc = 1:length(handles.sessionData{sessInd})
        % get the center points for each of the blocks in the original
        % image
        cfName = handles.nS.clusterInfoFileNameFcn(cc);
        clustInfoFile = load(cfName,'clusterBlockLocations');
        blockRanges = vertcat(clustInfoFile.clusterBlockLocations{:});
        blockRangeX = vertcat(blockRanges.blockRangeX);
        blockRangeY = vertcat(blockRanges.blockRangeY);
        % subtract offset
        y = (blockRangeY(:,1) + diff(blockRangeY,1,2)./2) - 1/2;
        x = (blockRangeX(:,1) + diff(blockRangeX,1,2)./2) - 1/2;
        
        
        pts = [x, y, squeeze( autoPts.xyzrcoClusterPeaks_auto(1:4,...
            chunkInd,:))'];
        % get the alignment points corresponding to those centerpoints
        handles.sessionData{sessInd}{chunkInd}.autoPts = num2cell(pts,2);
    end
end

% get the rois
if ~isfield(handles.sessionData{sessInd}{chunkInd},'roiMax') ...
        || isempty(handles.sessionData{sessInd}{chunkInd}.roiMax)
    
    try
        [handles.sessionData{sessInd}{chunkInd}.roiMax] = loadRoiMax(handles);
    catch
    end
    
end

set(handles.loadingText, 'Visible', 'off');
    
function [roiMaxToReturn] = loadRoiMax(handles)
%imageSize = ??handles.clustInfo.clusterParameters.imageSize

dateStr = handles.dateList{handles.currentSessionInd};
bsFileRoiMax = load(handles.nS.refLocBaselineFileName,'roiMax');
if ~isfield(bsFileRoiMax,'roiMax')
    bsFileRoiMax = {}; 
    roiMax = {}; 
else
    roiMax = bsFileRoiMax.roiMax;
end

for bb = 1:prod(handles.chunkParams{handles.currentSessionInd}{handles.currentChunkInd}.clusterParameters.nBlocks)
    
    % current using best angle as 0 even though we should really use the
    % true value
    bestAngle = 0; % REPLACE WITH REAL BEST ANGLE
    %bestAngle = handles.xyzrcoClusterPeaks(4,handles.currentChunkInd,bb);
    
    % if the rois are already in the localization file, just load them
    if ~isempty(bsFileRoiMax) && length(bsFileRoiMax.roiMax) >= handles.currentChunkInd && ...
            length(bsFileRoiMax.roiMax{handles.currentChunkInd}) >= bb 
        thisRoiMax = bsFileRoiMax.roiMax{handles.currentChunkInd}{bb};
        cImSz = size(thisRoiMax);
        %figure out how big cIm is going to be when rotated appropriately
        [rotH, rotW] = dimAfterRotation(cImSz(1), cImSz(2), ...
            bestAngle);
        roiMax{handles.currentChunkInd}{bb} = thisRoiMax;
    else
        % look for cell-finding file for this chunk
        cellfilename = sprintf(handles.nS.cellFileNameFcn(handles.currentChunkInd,bb),...
            dateStr);
        
        % load the rois from the cell file
        cellfile = load(cellfilename,'rois');

        thisRoiMax = max(bsxfun(@rdivide,cellfile.rois,max(max(cellfile.rois,[],1),[],2)),[],3);
        roiMax{handles.currentChunkInd}{bb} = thisRoiMax;
        save(handles.nS.refLocBaselineFileName, 'roiMax', '-append')
        
    end
    roiMaxToReturn = roiMax{handles.currentChunkInd};
end


function [roiMaxToReturn, roiMaxXRange, roiMaxYRange] = stitchRoisRigid(handles)
    bestAngle = 0;
    warning('best angle is hard coded as 0')

    sessInd = handles.currentSessionInd;
    chunkInd = handles.currentChunkInd;

    imageSize = handles.chunkParams{handles.currentSessionInd}{handles.currentChunkInd}.clusterParameters.imageSize;

    roiMaxXRange = [-handles.xPad imageSize(2)+handles.xPad];
    roiMaxYRange = [-handles.yPad imageSize(1)+handles.yPad];
    
    roiMaxToReturn = zeros(2.*[handles.yPad handles.xPad]+imageSize);

    
    for bb = 1:length(handles.sessionData{sessInd}{chunkInd}.roiMax)
        thisRoiMax = handles.sessionData{sessInd}{chunkInd}.roiMax{bb};
        cImSz = size(thisRoiMax);
        
        %figure out how big cIm is going to be when rotated appropriately
        [rotH, rotW] = dimAfterRotation(cImSz(1), cImSz(2), ...
            bestAngle);
        
        %pad roi weights with nans so they can rotate
        padSz = max([1 1], ceil([rotH-cImSz(1)  rotW-cImSz(2)]./2));
        
        % create the rotated image
        roiWPadded = padarray(thisRoiMax, padSz, nan,'both');
        roiWPadded = bsxfun(@rdivide,roiWPadded,max(max(roiWPadded,[],1),[],2));
        roiWLoc = true(size(roiWPadded,1), size(roiWPadded,2));
        roiWPaddedSz = size(roiWPadded);
        rotatedRoiW = rotateAndSelectBlock(max(roiWPadded,[],3),...
            roiWLoc,bestAngle);
        
        meanCenter = @(x)x-mean(x(:));
        
        % store in output variable
        autoPts = handles.sessionData{handles.currentSessionInd...
            }{handles.currentChunkInd}.autoPts;
        
        rotatedRoiX = ceil(repmat(autoPts{bb}(3) + meanCenter(1:roiWPaddedSz(2)),...
            roiWPaddedSz(1),1));
        rotatedRoiY = ceil(repmat(autoPts{bb}(4) + meanCenter(1:roiWPaddedSz(1))',...
            1, roiWPaddedSz(2)));
        rotatedRoiZ = repmat(autoPts{bb}(4), roiWPaddedSz(1), roiWPaddedSz(2));
        
        roiMaxToReturn(handles.yPad+(min(rotatedRoiY(:)):max(rotatedRoiY(:))),...
            handles.xPad+(min(rotatedRoiX(:)):max(rotatedRoiX(:)))) = ...
            (max(roiMaxToReturn(handles.yPad+(min(rotatedRoiY(:)):max(rotatedRoiY(:))),...
            handles.xPad+(min(rotatedRoiX(:)):max(rotatedRoiX(:)))),...
            rotatedRoiW));
    end

    
    
function [imToReturn, xRange, yRange] = stitchRoisNonRigid(handles)

tform = nonRigidTform(handles);
sessInd = handles.currentSessionInd;
chunkInd = handles.currentChunkInd;


for bb = 1:length(handles.sessionData{sessInd}{chunkInd}.roiMax)  
    thisBlockRois = handles.sessionData{sessInd}{chunkInd}.roiMax{bb};
    thisBlockXData = handles.sessionData{sessInd}{chunkInd}.autoPts{bb}(1) + ...
        [-1 1].*size(thisBlockRois,2)/2 + 1/2;
    thisBlockYData = handles.sessionData{sessInd}{chunkInd}.autoPts{bb}(2) + ...
        [-1 1].*size(thisBlockRois,1)/2 + 1/2;
    [imTransformed{bb}, xData, yData] = imtransform(thisBlockRois,...
        tform, 'udata',thisBlockXData,'vdata',thisBlockYData,'fillvalues', nan);
    
    % note the x and y values of pixels in the transformed image
    xVals{bb} = linspace(round(xData(1)),round(xData(2)),size(imTransformed{bb},2));
    yVals{bb} = linspace(round(yData(1)),round(yData(2)),size(imTransformed{bb},1));
    
    
end

xRange = [min([xVals{:}]) max([xVals{:}])];
yRange = [min([yVals{:}]) max([yVals{:}])];

for bb = 1:length(handles.sessionData{sessInd}{chunkInd}.roiMax)  
% set to nan where there is no image data
xInd = repmat(xVals{bb}-xRange(1)+1,length(yVals{bb}),1);
yInd = repmat(yVals{bb}-yRange(1)+1',1,length(xVals{bb}));
xInd(isnan(imTransformed{bb})) = nan;
yInd(isnan(imTransformed{bb})) = nan;

imToReturn(yInd, xInd) = imTransformed{bb};
end
                
                
                
                
function handles = makeBlockButtons(handles)
gridPos = get(handles.blockGridAx,'position');

handles.bg = uibuttongroup('Visible','off',...
                  'Position',gridPos);

              
% Create three radio buttons in the button group.
nX = handles.chunkParams{handles.currentSessionInd}{handles.currentChunkInd}.clusterParameters.nBlocks(1);
nY = handles.chunkParams{handles.currentSessionInd}{handles.currentChunkInd}.clusterParameters.nBlocks(2);
gapX = .1; gapY = .1; bounds = [0 0 1 1];
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
pause








% launch gui when click a block
function bselection(bg,event)
%     display(['Previous: ', (event.OldValue.String)]);
%     display(['Current: ', (event.NewValue.String)]);
%     display('------------------');
disp('hi')
% handles=guidata(bg);
% selectedBlock = str2num(event.Source.String)
% bs = load(sprintf(handles.nS.baselineNameFcn(handles.currentCluster,...
%     selectedBlock),...
%     handles.currentDateStr),'extraStats');
% cIm = normalizeToZeroOne(bs.extraStats.F_mean-handles.darkLevel);
% cImSz = size(cIm);%size(cIm);
% refIm = normalizeToZeroOne(double(handles.stack(:,:,handles.currentZ)));
% cpselect(cIm,refIm);
% % Update handles structure
% guidata(bg, handles);



% --- Executes just before reviewRegistration2 is made visible.
function reviewRegistration2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to reviewRegistration2 (see VARARGIN)


% variables that are currently hard-coded but shouldn't be
handles.yPad = 50;
handles.xPad = 50;

% Record mouse name
handles.mouse = varargin{1};

% Record recording site location
handles.location = varargin{2};

% optional stack
if length(varargin)>2
handles.reference = varargin{3};
end

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
updatePlot(handles)

% --- Executes on button press in rigidToggle.
function rigidToggle_Callback(hObject, eventdata, handles)
% hObject    handle to rigidToggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rigidToggle
updatePlot(handles)

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

% make new block buttons
%handles = makeBlockButtons(handles)

% turn hold off so that the figure won't keep it's zoom
% hold(handles.imRefAx,'off')


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

currProg = 1;
progStepSz = 20;
progSteps = linspace(1,length(handles.dateList),progStepSz);
fprintf('loading all datasets...')
for cc = 1:length(handles.nChunksInSession)
    handles.currentSessionInd = cc;
    
    % update the user on how the loading is going
    currProgNew = find(cc>=progSteps,1,'last');
    if currProgNew > currProg
        fprintf('%i%%...',currProg/progStepSz*100)
        currProg = currProgNew;
    end    
    
    % load each chunk separately
    for dd = 1:handles.nChunksInSession(cc)
        handles.currentChunkInd = dd;
        handles = loadDataset(handles);
    end
    
    guidata(hObject, handles)
    
end
sprintf('\n')



% --- Executes when selected object is changed in displayStyleGroup.
function displayStyleGroup_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in displayStyleGroup 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject,handles)
updatePlot(handles)
        

function tform = nonRigidTform(handles)

autoPts = cat(1,handles.sessionData{handles.currentSessionInd}{handles.currentChunkInd}.autoPts{:});
warning('doesn''t care about user selected points')

tform = cp2tform(autoPts(:,1:2), autoPts(:,3:4), 'lwm');


% 
% 
% 
% 
% % transform baseline image and each ROI (all at once)
% [dataTransformed, xData, yData] = imtransform(cat(3,baselineImage,cellfile.rois),tform,...
%     'udata',blockRange.blockRangeX,'vdata',blockRange.blockRangeY,'fillvalues', nan);
% 









function [] = updatePlot(handles)

cla(handles.imRefAx);

sessInd = handles.currentSessionInd;
chunkInd = handles.currentChunkInd;
    
autoPts = cat(1,handles.sessionData{sessInd}{chunkInd}.autoPts{:});

nBlocks = prod(handles.chunkParams{handles.currentSessionInd}{handles.currentChunkInd}.clusterParameters.nBlocks);

if get(handles.imgRefToggle,'Value') % show the image
    switch handles.displayStyleGroup.SelectedObject.String
        case 'Baseline'
            plot(handles.imRefAx,autoPts(:,3),autoPts(:,4),'.');
        case 'Blocks'
            plot(handles.imRefAx,autoPts(:,3),autoPts(:,4),'.');
        case 'ROIs'
            if get(handles.rigidToggle,'Value')
                % check to see if we need to create the rigid roi max
                if ~isfield(handles,'roiMaxRigid') | isempty(handles.roiMaxRigid{sessInd}{chunkInd})
                    [handles.sessionData{sessInd}{chunkInd}.roiMaxRigid,...
                        handles.sessionData{sessInd}{chunkInd}.roiMaxRigidXRange, ...
                        handles.sessionData{sessInd}{chunkInd}.roiMaxRigidYRange] = ...
                        stitchRoisRigid(handles);
                end
                % get the rigid roi max to plot it
                roiXData = handles.sessionData{sessInd}{chunkInd}.roiMaxRigidXRange(1):handles.sessionData{sessInd}{chunkInd}.roiMaxRigidXRange(end);
                roiYData = handles.sessionData{sessInd}{chunkInd}.roiMaxRigidYRange(1):handles.sessionData{sessInd}{chunkInd}.roiMaxRigidYRange(end);
                roisToPlot = handles.sessionData{sessInd}{chunkInd}.roiMaxRigid;
                
            else % use the nonrigid transformation
                % check to see if we need to create the nonrigid roi max
                if ~isfield(handles,'roiMaxNonRigid') | isempty(handles.roiMaxNonRigid{sessInd}{chunkInd})
                    [handles.sessionData{sessInd}{chunkInd}.roiMaxNonRigid,...
                        handles.sessionData{sessInd}{chunkInd}.roiMaxNonRigidXRange, ...
                        handles.sessionData{sessInd}{chunkInd}.roiMaxNonRigidYRange] = ...
                        stitchRoisNonRigid(handles);
                end
                % get the rigid roi max to plot it
                roiXData = handles.roiMaxNonRigidXRange{sessInd}{chunkInd}(1):handles.roiMaxNonRigidXRange{sessInd}{chunkInd}(end);
                roiYData = handles.roiMaxNonRigidYRange{sessInd}{chunkInd}(1):handles.roiMaxNonRigidYRange{sessInd}{chunkInd}(end);
                roisToPlot = handles.sessionData{sessInd}{chunkInd}.roiMaxNonRigid;
            end
            imagesc(roisToPlot, 'xdata', roiXData, 'ydata', roiYData, 'parent', handles.imRefAx);
            colormap(handles.imRefAx, gray)
            
    end
else
    modalZ = mode(autoPts(:,5));
    refXData = -handles.xPad:(handles.xPad+size(handles.reference,2));
    refYData = -handles.yPad:(handles.yPad+size(handles.reference,1));
    if get(handles.meanMultiZToggle, 'Value')
        refSlice = handles.reference(:,:,modalZ);
    else
        refSlice = handles.reference(:,:,modalZ);
    end
    paddedRef = padarray(refSlice,[handles.yPad handles.xPad],nan,'both');
    imagesc(paddedRef,'xdata',refXData,'ydata',refYData,'parent',handles.imRefAx);
    
end

hold(handles.imRefAx, 'on');


uicontrol(handles.sessionChunkSliderPage)
