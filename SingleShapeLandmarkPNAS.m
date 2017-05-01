function varargout = SingleShapeLandmarkPNAS(varargin)
% SINGLESHAPELANDMARKPNAS MATLAB code for SingleShapeLandmarkPNAS.fig
%      SINGLESHAPELANDMARKPNAS, by itself, creates a new SINGLESHAPELANDMARKPNAS or raises the existing
%      singleton*.
%
%      H = SINGLESHAPELANDMARKPNAS returns the handle to a new SINGLESHAPELANDMARKPNAS or the handle to
%      the existing singleton*.
%
%      SINGLESHAPELANDMARKPNAS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SINGLESHAPELANDMARKPNAS.M with the given input arguments.
%
%      SINGLESHAPELANDMARKPNAS('Property','Value',...) creates a new SINGLESHAPELANDMARKPNAS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SingleShapeLandmarkPNAS_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SingleShapeLandmarkPNAS_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SingleShapeLandmarkPNAS

% Last Modified: Tingran Gao (trgao10@math.duke.edu) April 13, 2017

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SingleShapeLandmarkPNAS_OpeningFcn, ...
                   'gui_OutputFcn',  @SingleShapeLandmarkPNAS_OutputFcn, ...
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


% --- Executes just before SingleShapeLandmarkPNAS is made visible.
function SingleShapeLandmarkPNAS_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SingleShapeLandmarkPNAS (see VARARGIN)

% Choose default command line output for SingleShapeLandmarkPNAS
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SingleShapeLandmarkPNAS wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SingleShapeLandmarkPNAS_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1

contents = cellstr(get(hObject,'String'));
load([hObject.UserData.sample_path contents{get(hObject,'Value')} '.mat']);

CameraUpVector = get(gca, 'CameraUpVector');
CameraPosition = get(gca, 'CameraPosition');
CameraTarget = get(gca, 'CameraTarget');
CameraViewAngle = get(gca, 'CameraViewAngle');

cla
color_data = repmat([0.9, 0.9, 0.8], G.nV, 1);
G.draw(struct('FaceColor', 'interp',...
    'FaceVertexCData', color_data, 'CDataMapping','scaled',...
    'EdgeColor', 'none', 'FaceAlpha', 1,...
    'AmbientStrength',0.3,'SpecularStrength',0.0));

set(gca, 'CameraUpVector',  CameraUpVector);
set(gca, 'CameraPosition',  CameraPosition);
set(gca, 'CameraTarget',    CameraTarget);
set(gca, 'CameraViewAngle', CameraViewAngle);

camlight('headlight');
camlight(180,0);


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

path(pathdef);
path(path, genpath('./utils'));

data_path = '~/Work/MATLAB/DATA/PNAS/';
sample_path = [data_path 'samples/'];
hObject.UserData.sample_path = sample_path;

FileList = getFileNames(sample_path);
set(hObject, 'String', FileList);

contents = cellstr(get(hObject,'String'));
load([sample_path contents{get(hObject,'Value')} '.mat']);

color_data = repmat([0.9, 0.9, 0.8], G.nV, 1);
G.draw(struct('FaceColor', 'interp',...
    'FaceVertexCData', color_data, 'CDataMapping','scaled',...
    'EdgeColor', 'none', 'FaceAlpha', 1,...
    'AmbientStrength',0.3,'SpecularStrength',0.0));

camlight('headlight');
camlight(180,0);


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

userdata = get(gcf,'UserData');
G = userdata.mesh;
if ~isfield(G.Aux,'fullPhi')
    [Cgauss,Cmean,Umin,Umax,Cmin,Cmax,Normal] = G.ComputeCurvature();
    DNE = Cmin.^2+Cmax.^2;
    DNE(DNE<median(DNE)) = 0;
    [EdgeIdxI,EdgeIdxJ] = find(tril(G.E));
    bandwidth = mean(sqrt(sum((G.V(:,EdgeIdxI)-G.V(:,EdgeIdxJ)).^2)))/4;
    PDistMat = squareform(pdist(G.V'));
    fullPhi = exp(-PDistMat.^2/bandwidth);
%     VertAreaMeas = diag(G.Aux.VertArea);
%     VertAreaMeas = diag(G.Aux.VertArea)+diag(abs(Cgauss)/sum(abs(Cgauss))+abs(Cmean)/sum(abs(Cmean))+sqrt(DNE)/sum(sqrt(DNE)));
    VertAreaMeas = diag(G.Aux.VertArea)+diag(abs(Cgauss)/sum(abs(Cgauss))+abs(Cmean)/sum(abs(Cmean)));
    userdata.mesh.Aux.fullPhi = fullPhi;
    userdata.mesh.Aux.VertAreaMeas = VertAreaMeas;
else
    fullPhi = userdata.mesh.Aux.fullPhi;
    VertAreaMeas = userdata.mesh.Aux.VertAreaMeas;
end

if ~isfield(G.Aux, 'GPLmkIdx')
    ptuq = sum(fullPhi.*(VertAreaMeas*fullPhi))';
%     InitLmkIdx = G.Aux.ConfMaxInds;
    [~,InitLmkIdx] = max(ptuq);
    userdata.mesh.Aux.GPLmkIdx = InitLmkIdx;
    
    CameraUpVector = get(gca, 'CameraUpVector');
    CameraPosition = get(gca, 'CameraPosition');
    CameraTarget = get(gca, 'CameraTarget');
    CameraViewAngle = get(gca, 'CameraViewAngle');
    
    cla
    G.ViewFunctionOnMesh(ptuq,struct('mode','native'));
    hold on
    scatter3(G.V(1,InitLmkIdx),G.V(2,InitLmkIdx),G.V(3,InitLmkIdx),20,'r','filled');
    
    set(gca, 'CameraUpVector',  CameraUpVector);
    set(gca, 'CameraPosition',  CameraPosition);
    set(gca, 'CameraTarget',    CameraTarget);
    set(gca, 'CameraViewAngle', CameraViewAngle);
else
    lmkPhi = fullPhi(:,userdata.mesh.Aux.GPLmkIdx);
    projNullPhi = VertAreaMeas-VertAreaMeas*lmkPhi*((lmkPhi'*VertAreaMeas*lmkPhi)\(lmkPhi'*VertAreaMeas));
    ptuq = sum(fullPhi.*(projNullPhi*fullPhi))';
    [~,maxUQIdx] = max(ptuq);
    userdata.mesh.Aux.GPLmkIdx = [userdata.mesh.Aux.GPLmkIdx;maxUQIdx];
    
    CameraUpVector = get(gca, 'CameraUpVector');
    CameraPosition = get(gca, 'CameraPosition');
    CameraTarget = get(gca, 'CameraTarget');
    CameraViewAngle = get(gca, 'CameraViewAngle');
    
    cla
    G.ViewFunctionOnMesh(ptuq,struct('mode','native'));
    hold on
    scatter3(G.V(1,userdata.mesh.Aux.GPLmkIdx(1:end-1)),...
        G.V(2,userdata.mesh.Aux.GPLmkIdx(1:end-1)),...
        G.V(3,userdata.mesh.Aux.GPLmkIdx(1:end-1)),20,'r','filled');
    scatter3(G.V(1,userdata.mesh.Aux.GPLmkIdx(end)),...
        G.V(2,userdata.mesh.Aux.GPLmkIdx(end)),...
        G.V(3,userdata.mesh.Aux.GPLmkIdx(end)),20,'g','filled');
    
    set(gca, 'CameraUpVector',  CameraUpVector);
    set(gca, 'CameraPosition',  CameraPosition);
    set(gca, 'CameraTarget',    CameraTarget);
    set(gca, 'CameraViewAngle', CameraViewAngle);
end


