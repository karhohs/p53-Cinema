function varargout = manualTracking(varargin)
% MANUALTRACKING MATLAB code for manualTracking.fig
%      MANUALTRACKING, by itself, creates a new MANUALTRACKING or raises the existing
%      singleton*.
%
%      H = MANUALTRACKING returns the handle to a new MANUALTRACKING or the handle to
%      the existing singleton*.
%
%      MANUALTRACKING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MANUALTRACKING.M with the given input arguments.
%
%      MANUALTRACKING('Property','Value',...) creates a new MANUALTRACKING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before manualTracking_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to manualTracking_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help manualTracking

% Last Modified by GUIDE v2.5 21-Sep-2012 16:15:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @manualTracking_OpeningFcn, ...
                   'gui_OutputFcn',  @manualTracking_OutputFcn, ...
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


% --- Executes just before manualTracking is made visible.
function manualTracking_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to manualTracking (see VARARGIN)

% Choose default command line output for manualTracking
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

handles.ind = 1;
guidata(hObject, handles);
% UIWAIT makes manualTracking wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = manualTracking_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function editLogPath_Callback(hObject, eventdata, handles) %#ok<*INUSD>
% hObject    handle to editLogPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editLogPath as text
%        str2double(get(hObject,'String')) returns contents of editLogPath as a double


% --- Executes during object creation, after setting all properties.
function editLogPath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editLogPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonLogPath.
function pushbuttonLogPath_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonLogPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

folder_name = uigetdir;
set(handles.editLogPath,'String',folder_name);



function editStackPath_Callback(hObject, eventdata, handles)
% hObject    handle to editStackPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editStackPath as text
%        str2double(get(hObject,'String')) returns contents of editStackPath as a double


% --- Executes during object creation, after setting all properties.
function editStackPath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editStackPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonStackPath.
function pushbuttonStackPath_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonStackPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

folder_name = uigetdir;
set(handles.editStackPath,'String',folder_name);


% --- Executes on button press in pushbuttonExtractData.
function pushbuttonExtractData_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonExtractData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.textStepTwoFinished, 'Visible', 'off');
pause(0.1); %this 
logpath = get(handles.editLogPath, 'String');
stackpath = get(handles.editStackPath, 'String');
fluorchan = get(handles.editFluorescentChannels, 'String'); %assume CSV
C = textscan(fluorchan, '%s', 'delimiter', ', ', 'MultipleDelimsAsOne', 1);
for i = 1:length(C{1})
    processManualSegTrackViaImageJ(logpath, stackpath, 'fluorchan', C{1}{i},'phaseratio',4)
end
c = clock;
set(handles.textStepTwoFinished, 'String', sprintf('Finished! @ %02d:%02d',c(4),c(5)));
set(handles.textStepTwoFinished, 'Visible', 'on');
str = sprintf('dynamics%s', C{1}{1});
set(handles.editDataPath, 'String', fullfile(logpath,str));




function editFluorescentChannels_Callback(hObject, eventdata, handles)
% hObject    handle to editFluorescentChannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFluorescentChannels as text
%        str2double(get(hObject,'String')) returns contents of editFluorescentChannels as a double


% --- Executes during object creation, after setting all properties.
function editFluorescentChannels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFluorescentChannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonPrev.
function pushbuttonPrev_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonPrev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.ind~=1
   handles.ind = handles.ind - 1;
end
set(handles.editCurrentCell,'String',num2str(handles.ind));
plot(handles.axes1, handles.data(handles.ind).meanIntensity);
guidata(hObject, handles);

% --- Executes on button press in pushbuttonNext.
function pushbuttonNext_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonNext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.ind~=length(handles.data)
   handles.ind = handles.ind + 1;
end
set(handles.editCurrentCell,'String',num2str(handles.ind));
plot(handles.axes1, handles.data(handles.ind).meanIntensity);
guidata(hObject, handles);


function editCurrentCell_Callback(hObject, eventdata, handles)
% hObject    handle to editCurrentCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editCurrentCell as text
%        str2double(get(hObject,'String')) returns contents of editCurrentCell as a double
currentInd = str2double(get(handles.editCurrentCell,'String'));
handles.ind = currentInd;
plot(handles.axes1, handles.data(handles.ind).meanIntensity);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function editCurrentCell_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCurrentCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editDataPath_Callback(hObject, eventdata, handles)
% hObject    handle to editDataPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDataPath as text
%        str2double(get(hObject,'String')) returns contents of editDataPath as a double


% --- Executes during object creation, after setting all properties.
function editDataPath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDataPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonDataPath.
function pushbuttonDataPath_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonDataPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, foldername, ~] = uigetfile;
set(handles.editDataPath,'String',fullfile(foldername, filename));

% --- Executes on button press in pushbuttonLoadData.
function pushbuttonLoadData_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonLoadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
datapath = get(handles.editDataPath, 'String');
load(datapath);
c = clock;
set(handles.textLoadData, 'String', sprintf('Data Loaded! @ %02d:%02d',c(4),c(5)));
set(handles.textLoadData, 'Visible', 'on');
handles.data = unitOfLife;
guidata(hObject, handles);



function editRatio_Callback(hObject, eventdata, handles)
% hObject    handle to editRatio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editRatio as text
%        str2double(get(hObject,'String')) returns contents of editRatio as a double


% --- Executes during object creation, after setting all properties.
function editRatio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editRatio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
