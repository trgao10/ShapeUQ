function varargout = InteractiveSingleShapeLandmarkPNAS(varargin)
% INTERACTIVESINGLESHAPELANDMARKPNAS MATLAB code for InteractiveSingleShapeLandmarkPNAS.fig
%      INTERACTIVESINGLESHAPELANDMARKPNAS, by itself, creates a new INTERACTIVESINGLESHAPELANDMARKPNAS or raises the existing
%      singleton*.
%
%      H = INTERACTIVESINGLESHAPELANDMARKPNAS returns the handle to a new INTERACTIVESINGLESHAPELANDMARKPNAS or the handle to
%      the existing singleton*.
%
%      INTERACTIVESINGLESHAPELANDMARKPNAS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INTERACTIVESINGLESHAPELANDMARKPNAS.M with the given input arguments.
%
%      INTERACTIVESINGLESHAPELANDMARKPNAS('Property','Value',...) creates a new INTERACTIVESINGLESHAPELANDMARKPNAS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before InteractiveSingleShapeLandmarkPNAS_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to InteractiveSingleShapeLandmarkPNAS_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help InteractiveSingleShapeLandmarkPNAS

% Last Modified: Tingran Gao (trgao10@math.duke.edu) April 13, 2017

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @InteractiveSingleShapeLandmarkPNAS_OpeningFcn, ...
                   'gui_OutputFcn',  @InteractiveSingleShapeLandmarkPNAS_OutputFcn, ...
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


% --- Executes just before InteractiveSingleShapeLandmarkPNAS is made visible.
function InteractiveSingleShapeLandmarkPNAS_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to InteractiveSingleShapeLandmarkPNAS (see VARARGIN)

% Choose default command line output for InteractiveSingleShapeLandmarkPNAS
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes InteractiveSingleShapeLandmarkPNAS wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = InteractiveSingleShapeLandmarkPNAS_OutputFcn(hObject, eventdata, handles) 
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

userdata = get(gcf,'UserData');
userdata.ToggleInteractiveLandmarking = 'off';
set(gcf, 'userdata', userdata);
set(gcf, 'KeyPressFcn', {@ToggleInteractiveLandmarking});
set(gcf, 'Name', 'InteractiveSingleShapeLandmarkPNAS');


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

userdata = get(gcf,'UserData');
userdata.ToggleInteractiveLandmarking = 'off';
userdata.currentWindowButtonDownFcn = get(gcf, 'WindowButtonDownFcn');
set(gcf, 'userdata', userdata);
set(gcf, 'KeyPressFcn', {@ToggleInteractiveLandmarking});
set(gcf, 'Name', 'InteractiveSingleShapeLandmarkPNAS');


function ToggleInteractiveLandmarking(hObject, eventData)
    
if(strcmp(eventData.Key, 'return'))
    userdata = get(gcf, 'userdata');
    if strcmpi(userdata.ToggleInteractiveLandmarking,'off')
        userdata.currentWindowButtonDownFcn = get(gcf, 'WindowButtonDownFcn');
        userdata.ToggleInteractiveLandmarking = 'on';
        frameTitleStr = get(gcf, 'Name');
        if isempty(strfind(lower(frameTitleStr), 'off')) && isempty(strfind(lower(frameTitleStr), 'on'))
            frameTitleStr = [frameTitleStr ' (Landmarking Mode - ON)'];
            
            G = userdata.mesh;
            if ~isfield(G.Aux,'fullPhi')
                [Cgauss,Cmean,Umin,Umax,Cmin,Cmax,Normal] = G.ComputeCurvature();
                DNE = Cmin.^2+Cmax.^2;
                DNE(DNE<median(DNE)) = 0;
                [EdgeIdxI,EdgeIdxJ] = find(tril(G.E));
                bandwidth = mean(sqrt(sum((G.V(:,EdgeIdxI)-G.V(:,EdgeIdxJ)).^2)))/4;
                PDistMat = squareform(pdist(G.V'));
                fullPhi = exp(-PDistMat.^2/bandwidth);
                VertAreaMeas = diag(G.Aux.VertArea)+diag(abs(Cgauss)/sum(abs(Cgauss))+abs(Cmean)/sum(abs(Cmean)));
                userdata.mesh.Aux.fullPhi = fullPhi;
                userdata.mesh.Aux.VertAreaMeas = VertAreaMeas;
            else
                fullPhi = userdata.mesh.Aux.fullPhi;
                VertAreaMeas = userdata.mesh.Aux.VertAreaMeas;
            end
            
            if ~isfield(G.Aux, 'GPLmkIdx')
                ptuq = sum(fullPhi.*(VertAreaMeas*fullPhi))';
            else
                lmkPhi = fullPhi(:,userdata.mesh.Aux.GPLmkIdx);
                projNullPhi = VertAreaMeas-VertAreaMeas*lmkPhi*((lmkPhi'*VertAreaMeas*lmkPhi)\(lmkPhi'*VertAreaMeas));
                ptuq = sum(fullPhi.*(projNullPhi*fullPhi))';
            end
            
            CameraUpVector = get(gca, 'CameraUpVector');
            CameraPosition = get(gca, 'CameraPosition');
            CameraTarget = get(gca, 'CameraTarget');
            CameraViewAngle = get(gca, 'CameraViewAngle');
            
            cla
            G.ViewFunctionOnMesh(ptuq,struct('mode','native'));
            set(gcf, 'KeyPressFcn', {@ToggleInteractiveLandmarking});
            if isfield(G.Aux, 'GPLmkIdx')
                if length(G.Aux.GPLmkIdx) == 1
                    scatter3(G.V(1,G.Aux.GPLmkIdx),G.V(2,G.Aux.GPLmkIdx),G.V(3,G.Aux.GPLmkIdx),20,'g','filled');
                else
                    scatter3(G.V(1,userdata.mesh.Aux.GPLmkIdx(1:end-1)),...
                        G.V(2,userdata.mesh.Aux.GPLmkIdx(1:end-1)),...
                        G.V(3,userdata.mesh.Aux.GPLmkIdx(1:end-1)),20,'r','filled');
                    scatter3(G.V(1,userdata.mesh.Aux.GPLmkIdx(end)),...
                        G.V(2,userdata.mesh.Aux.GPLmkIdx(end)),...
                        G.V(3,userdata.mesh.Aux.GPLmkIdx(end)),20,'g','filled');
                end
            end
            
            set(gca, 'CameraUpVector',  CameraUpVector);
            set(gca, 'CameraPosition',  CameraPosition);
            set(gca, 'CameraTarget',    CameraTarget);
            set(gca, 'CameraViewAngle', CameraViewAngle);
        else
           frameTitleStr = strrep(frameTitleStr, 'OFF', 'ON');
        end
        
        set(gcf, 'WindowButtonDownFcn', {@tgDataCursor});
        set(gcf, 'Name', frameTitleStr);
        set(gcf, 'userdata', userdata);
    elseif (strcmpi(userdata.ToggleInteractiveLandmarking,'on'))
        set(gcf, 'WindowButtonDownFcn', userdata.currentWindowButtonDownFcn);
        userdata.ToggleInteractiveLandmarking = 'off';
        frameTitleStr = get(gcf, 'Name');
        frameTitleStr = strrep(frameTitleStr, 'ON', 'OFF');
        set(gcf, 'Name', frameTitleStr);
        set(gcf, 'userdata', userdata);
    else
        disp('Something is not right...');
    end
end


function tgDataCursor(hObject, eventData)

userdata = get(gcf, 'userdata');
pt = userdata.mesh.LineMeshIntersect(get(gca, 'CurrentPoint'));

if ~isfield(userdata.mesh.Aux, 'kdtree')
    userdata.mesh.Aux.kdtree = kdtree_build(userdata.mesh.V');
    set(gcf, 'userdata', userdata);
end

if ~isempty(pt)
    
    nearestVertex = kdtree_nearest_neighbor(userdata.mesh.Aux.kdtree,pt');
    InitLmkIdx = nearestVertex;
    
    disp(nearestVertex);
    
    G = userdata.mesh;
    fullPhi = userdata.mesh.Aux.fullPhi;
    VertAreaMeas = userdata.mesh.Aux.VertAreaMeas;
    
    if ~isfield(G.Aux, 'GPLmkIdx')
        userdata.mesh.Aux.GPLmkIdx = InitLmkIdx;
        
        lmkPhi = fullPhi(:,userdata.mesh.Aux.GPLmkIdx);
        projNullPhi = VertAreaMeas-VertAreaMeas*lmkPhi*((lmkPhi'*VertAreaMeas*lmkPhi)\(lmkPhi'*VertAreaMeas));
        ptuq = sum(fullPhi.*(projNullPhi*fullPhi))';
        
        CameraUpVector = get(gca, 'CameraUpVector');
        CameraPosition = get(gca, 'CameraPosition');
        CameraTarget = get(gca, 'CameraTarget');
        CameraViewAngle = get(gca, 'CameraViewAngle');
        
        cla
        G.ViewFunctionOnMesh(ptuq,struct('mode','native'));
        hold on
        scatter3(G.V(1,InitLmkIdx),G.V(2,InitLmkIdx),G.V(3,InitLmkIdx),20,'g','filled');
        
        set(gcf, 'KeyPressFcn', {@ToggleInteractiveLandmarking});
        set(gcf, 'WindowButtonDownFcn', {@tgDataCursor});
        
        set(gca, 'CameraUpVector',  CameraUpVector);
        set(gca, 'CameraPosition',  CameraPosition);
        set(gca, 'CameraTarget',    CameraTarget);
        set(gca, 'CameraViewAngle', CameraViewAngle);
    else
        userdata.mesh.Aux.GPLmkIdx = [userdata.mesh.Aux.GPLmkIdx;InitLmkIdx];
        
        lmkPhi = fullPhi(:,userdata.mesh.Aux.GPLmkIdx);
        projNullPhi = VertAreaMeas-VertAreaMeas*lmkPhi*((lmkPhi'*VertAreaMeas*lmkPhi)\(lmkPhi'*VertAreaMeas));
        ptuq = sum(fullPhi.*(projNullPhi*fullPhi))';
        
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
        
        set(gcf, 'KeyPressFcn', {@ToggleInteractiveLandmarking});
        set(gcf, 'WindowButtonDownFcn', {@tgDataCursor});
        
        set(gca, 'CameraUpVector',  CameraUpVector);
        set(gca, 'CameraPosition',  CameraPosition);
        set(gca, 'CameraTarget',    CameraTarget);
        set(gca, 'CameraViewAngle', CameraViewAngle);
    end

end

set(gcf, 'userdata', userdata);
set(gcf, 'KeyPressFcn', {@ToggleInteractiveLandmarking});


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

userdata.ToggleInteractiveLandmarking = 'off';
set(gcf, 'userdata', userdata);
set(gcf, 'Name', 'InteractiveSingleShapeLandmarkPNAS');
set(gcf, 'KeyPressFcn', {@ToggleInteractiveLandmarking});


% --- Executes on key press with focus on pushbutton1 and none of its controls.
function pushbutton1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

userdata = get(gcf,'UserData');
if strcmpi(userdata.ToggleInteractiveLandmarking, 'off')
    userdata.ToggleInteractiveLandmarking = 'on';
    set(gcf, 'Name', 'InteractiveSingleShapeLandmarkPNAS (Landmarking Mode - ON)');
    set(gcf, 'KeyPressFcn', {@ToggleInteractiveLandmarking});
    set(gcf, 'WindowButtonDownFcn', {@tgDataCursor});    
else
    userdata.ToggleInteractiveLandmarking = 'off';
    set(gcf, 'Name', 'InteractiveSingleShapeLandmarkPNAS (Landmarking Mode - OFF)');
    set(gcf, 'KeyPressFcn', {@ToggleInteractiveLandmarking});
    set(gcf, 'WindowButtonDownFcn', userdata.currentWindowButtonDownFcn);
end

set(gcf, 'userdata', userdata);
