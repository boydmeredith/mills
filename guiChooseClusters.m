function varargout = guiChooseClusters(varargin)
% GUICHOOSECLUSTERS MATLAB code for guiChooseClusters.fig
%      GUICHOOSECLUSTERS, by itself, creates a new GUICHOOSECLUSTERS or raises the existing
%      singleton*.
%
%      H = GUICHOOSECLUSTERS returns the handle to a new GUICHOOSECLUSTERS or the handle to
%      the existing singleton*.
%
%      GUICHOOSECLUSTERS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUICHOOSECLUSTERS.M with the given input arguments.
%
%      GUICHOOSECLUSTERS('Property','Value',...) creates a new GUICHOOSECLUSTERS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before guiChooseClusters_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to guiChooseClusters_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help guiChooseClusters

% Last Modified by GUIDE v2.5 03-Aug-2016 19:31:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @guiChooseClusters_OpeningFcn, ...
                   'gui_OutputFcn',  @guiChooseClusters_OutputFcn, ...
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


% --- Executes just before guiChooseClusters is made visible.
function guiChooseClusters_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to guiChooseClusters (see VARARGIN)

% Choose default command line output for guiChooseClusters
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes guiChooseClusters wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = guiChooseClusters_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in buttonReset.
function buttonReset_Callback(hObject, eventdata, handles)
% hObject    handle to buttonReset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

plot(handles.axesMain,rand(100,3))