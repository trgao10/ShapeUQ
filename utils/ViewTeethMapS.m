function ViewTeethMapS(GM,GN,maps,options)
%VIEWTEMPLATES: visualize a pair of maps between GM and GN
%   GM, GN:     Mesh objects
%   maps:       1x2 cell array containing map12 (GM->GN) and map21 (GN->GM)
%   options:
%      Type:            {'full'} | 'sample'
%      Landmarks:       {'on'}  | 'off'
%      LandmarksPath:   path to landmarks_teeth.mat
%      MeshesPath:      path to the folder containing all meshes
%
%   Tingran Gao, trgao10@math.duke.edu
%   last modified: 15 Aug 2014
%

if nargin<4
    options = [];
end
Type = getoptions(options,'Type','full');
Landmarks = getoptions(options,'Landmarks','on');

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% teeth panel
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if strcmpi(Landmarks, 'on')
    figure('Unit', 'pixel', 'Position', [0,0,1200,800], 'Tag', 'TEETH');
    NumRows = 2;
    NumCols = 2;
    h = zeros(1,4);
else
    figure('Unit', 'pixel', 'Position', [0,0,1200,400], 'Tag', 'TEETH');
    NumRows = 1;
    NumCols = 2;
    h = zeros(1,2);
end
set(gcf, 'ToolBar', 'none');

mesh_list = cell(1,2);

mesh_list{1} = GM;
h(1) = subplot(NumRows,NumCols,1);
GM.draw(struct('FaceColor',[0.9 0.9 0.8],'EdgeColor','none','FaceAlpha',1,'AmbientStrength',0.3,'SpecularStrength',0.0));
set(gca,'Tag','Source');
hold on;
camlight('headlight');
camlight(180,0);
lighting phong;
title(GM.Aux.name);
if (strcmpi(Type, 'sample'))
    samples = GM.V(:,GM.Aux.VertSampInd);
    scatter3(samples(1,:),samples(2,:),samples(3,:),5,'k','filled');
end
if (strcmpi(Landmarks, 'on'))
    % extract landmarks
    [IndsOnSource, ~] = GetLandmarks(GM.Aux.name,options.LandmarksPath,[options.MeshesPath GM.Aux.name options.MeshSuffix],options);
    draw_landmarks(GM.V,IndsOnSource);
end

mesh_list{2} = GN;
h(2) = subplot(NumRows,NumCols,2);
GN.draw(struct('FaceColor',[0.9 0.9 0.8],'EdgeColor','none','FaceAlpha',1,'AmbientStrength',0.3,'SpecularStrength',0.0));
set(gca,'Tag','Target');
hold on;
camlight('headlight');
camlight(180,0);
lighting phong;
title(GN.Aux.name);
if (strcmpi(Type, 'sample'))
    samples = GN.V(:,GN.Aux.VertSampInd);
    scatter3(samples(1,:),samples(2,:),samples(3,:),5,'k','filled');
end
if (strcmpi(Landmarks, 'on'))
    % extract landmarks
    [IndsOnTarget, ~] = GetLandmarks(GN.Aux.name,options.LandmarksPath,[options.MeshesPath GN.Aux.name options.MeshSuffix],options);
    draw_landmarks(GN.V,IndsOnTarget);
end

if (strcmpi(Landmarks, 'on'))
    h(3) = subplot(NumRows,NumCols,3);
    GM.draw(struct('FaceColor',[0.9 0.9 0.8],'EdgeColor','none','FaceAlpha',1,'AmbientStrength',0.3,'SpecularStrength',0.0));
    set(gca,'Tag','PropagSource');
    hold on;
    camlight('headlight');
    camlight(180,0);
    lighting phong;
    if (strcmpi(Type, 'sample'))
        samples = GM.V(:,GM.Aux.VertSampInd);
        scatter3(samples(1,:),samples(2,:),samples(3,:),5,'k','filled');
    end
    map = maps{2};
    PropagLandmarkCoords = GM.V(:,map(IndsOnTarget));
    draw_landmarks(GM.V,PropagLandmarkCoords,'Coords');

    h(4) = subplot(NumRows,NumCols,4);
    GN.draw(struct('FaceColor',[0.9 0.9 0.8],'EdgeColor','none','FaceAlpha',1,'AmbientStrength',0.3,'SpecularStrength',0.0));
    set(gca,'Tag','PropagTarget');
    hold on;
    camlight('headlight');
    camlight(180,0);
    lighting phong;
    if (strcmpi(Type, 'sample'))
        samples = GM.V(:,GM.Aux.VertSampInd);
        scatter3(samples(1,:),samples(2,:),samples(3,:),5,'k','filled');
    end
    map = maps{1};
    PropagLandmarkCoords = GN.V(:,map(IndsOnSource));
    draw_landmarks(GN.V,PropagLandmarkCoords,'Coords');
end

Link = linkprop(h, {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'CameraViewAngle'});
setappdata(gcf, 'StoreTheLink', Link);

set(gca, 'CameraUpVector', [0.8469,-0.5272,-0.0696]);
set(gca, 'CameraPosition', [0.0425,0.5383,-3.7461]);
set(gca, 'CameraTarget', [0.0122,-0.0075,0.0173]);
set(gca, 'CameraViewAngle', 10.5477);

userdata.maps = maps;
userdata.Type = Type;
userdata.mesh_list = mesh_list;
userdata.ToggleOnTeethSelect3DFlag = 'off';
set(gcf, 'userdata', userdata);
set(gcf, 'KeyPressFcn', {@ToggleOnTeethSelect3D});

end

% function [Inds, Coords] = ExtractLandmarks(GM, options)
% 
% Type = getoptions(options,'Type','full');
% 
% LandmarkFile = load(options.LandmarksPath);
% rawLandmarks = LandmarkFile.PP(strcmpi(LandmarkFile.names, GM.Aux.name),1:options.NumLandmark,:);
% Landmarks = zeros(size(rawLandmarks,2),3);
% for k=1:size(rawLandmarks,2)
%     Landmarks(k,:) = [rawLandmarks(1,k,1), rawLandmarks(1,k,2), rawLandmarks(1,k,3)];
% end
% [V,F] = read_off([options.MeshesPath GM.Aux.name options.MeshSuffix]);
% V = V';
% F = F';
% area = CORR_calculate_area(F,V);
% V = V*sqrt(1/area);
% tree = KDTreeSearcher(V);
% LandmarkVertInds = tree.knnsearch(Landmarks);
% 
% if strcmpi(Type, 'sample')
%     LandmarkOnSampleInds = GM.Aux.V2S(LandmarkVertInds);
%     LandmarkOnSamples = GM.V(:,LandmarkOnSampleInds);
%     
%     LandmarkToNearestSample = zeros(size(LandmarkOnSampleInds));
%     for k=1:length(LandmarkToNearestSample)
%         LandmarkToNearestSample(k) = find(GM.Aux.VertSampInd==LandmarkOnSampleInds(k));
%     end
%     Inds = LandmarkToNearestSample;
%     Coords = LandmarkOnSamples;
% elseif strcmpi(Type, 'full')
%     Inds = LandmarkVertInds;
%     Coords = GM.V(:,Inds);
% end
%
% end

function draw_landmarks(V,IndsOnSource,Type)

if nargin<3
    Type = 'Inds';
end

if strcmpi(Type,'Coords')
    Landmarks = IndsOnSource;
else
    if (size(V,1)==3)
        V = V';
    end
    Landmarks = V(IndsOnSource,:)';
end

NumLandmarks = max(size(Landmarks));
    
colmap =  [1,0,0;0,1,0;0,0,1;1,1,0;1,0,1;0,1,1];
colmap = [colmap;colmap*0.5;colmap*0.25];

R = 0.015;
CORR_draw_spheres(Landmarks',R,colmap(1:NumLandmarks,:));

end

function ToggleOnTeethSelect3D(hObject, eventData)

if(strcmp(eventData.Key, 'space'))
    userdata = get(gcf, 'userdata');
    if strcmpi(userdata.ToggleOnTeethSelect3DFlag,'off')
        userdata.currentWindowButtonDownFcn = get(gcf, 'WindowButtonDownFcn');
        set(gcf, 'WindowButtonDownFcn',...
            {@callbackOnTeethDataCursor, userdata.mesh_list, userdata.Type, userdata.maps});
        userdata.ToggleOnTeethSelect3DFlag = 'on';
        set(gcf, 'Name',...
            ['OnTeethSelect3D ' userdata.ToggleOnTeethSelect3DFlag]);
        set(gcf, 'userdata', userdata);
    elseif (strcmpi(userdata.ToggleOnTeethSelect3DFlag,'on'))
        set(gcf, 'WindowButtonDownFcn',...
            userdata.currentWindowButtonDownFcn);
        userdata.ToggleOnTeethSelect3DFlag = 'off';
        set(gcf, 'Name',...
            ['OnTeethSelect3D ' userdata.ToggleOnTeethSelect3DFlag]);
        set(gcf, 'userdata', userdata);
    else
        disp('Something is not right...');
    end
end

end

function callbackOnTeethDataCursor(hObject,eventData,mesh_list,Type,maps)

MarkerSize = 30;

[~, selectedPoint, ~] = select3d;

if (isempty(selectedPoint))
    return;
end

h = findall(gcf,'Tag','pt');
delete(h);

if strcmpi(get(gca,'Tag'),'Source')||strcmpi(get(gca,'Tag'),'PropagSource')
    map = maps{1};
    
    GM = mesh_list{1};
    if (strcmpi(Type, 'sample'))
        tree = KDTreeSearcher(GM.V(:,GM.Aux.VertSampInd)');
        pointInSourceSampleIndex = tree.knnsearch(selectedPoint');
        pointOnSourceIndex = GM.Aux.VertSampInd(pointInSourceSampleIndex);
        pointOnSource = GM.V(:,pointOnSourceIndex);
    elseif (strcmpi(Type, 'full'))
        tree = KDTreeSearcher(GM.V');
        pointOnSourceIndex = tree.knnsearch(selectedPoint');
        pointOnSource = GM.V(:,pointOnSourceIndex);
    end
    clear GM;
    h = plot3(findobj('Tag','Source'), pointOnSource(1,:),...
        pointOnSource(2,:), pointOnSource(3,:), 'r.', 'MarkerSize', MarkerSize);
    set(h,'Tag','pt');
    if ~isempty(findobj('Tag','PropagSource'))
        h = plot3(findobj('Tag','PropagSource'), pointOnSource(1,:),...
            pointOnSource(2,:), pointOnSource(3,:), 'r.', 'MarkerSize', MarkerSize);
        set(h,'Tag','pt');
    end
    
    GM = mesh_list{2};
    if (strcmpi(Type, 'sample'))
        tree = KDTreeSearcher(GM.V(:,GM.Aux.VertSampInd)');
        pointOnTargetIndex = map(pointOnSourceIndex);
        pointOnTarget = GM.V(:,pointOnTargetIndex);
        pointOnTargetSampleIndex = tree.knnsearch(pointOnTarget');
        pointOnTargetIndex = GM.Aux.VertSampInd(pointOnTargetSampleIndex);
        pointOnTarget = GM.V(:,pointOnTargetIndex);
    elseif (strcmpi(Type, 'full'))
        pointOnTarget = GM.V(:,map(pointOnSourceIndex));
    end
    clear GM;
    h = plot3(findobj('Tag','Target'), pointOnTarget(1,:),...
        pointOnTarget(2,:), pointOnTarget(3,:), 'g.', 'MarkerSize', MarkerSize);
    set(h,'Tag','pt');
    if ~isempty(findobj('Tag','PropagTarget'))
        h = plot3(findobj('Tag','PropagTarget'), pointOnTarget(1,:),...
            pointOnTarget(2,:), pointOnTarget(3,:), 'g.', 'MarkerSize', MarkerSize);
        set(h,'Tag','pt');
    end
    
elseif strcmpi(get(gca,'Tag'),'Target')||strcmpi(get(gca,'Tag'),'PropagTarget')
    map = maps{2};
    
    GM = mesh_list{2};
    if (strcmpi(Type, 'sample'))
        tree = KDTreeSearcher(GM.V(:,GM.Aux.VertSampInd)');
        pointInSourceSampleIndex = tree.knnsearch(selectedPoint');
        pointOnSourceIndex = GM.Aux.VertSampInd(pointInSourceSampleIndex);
        pointOnSource = GM.V(:,pointOnSourceIndex);
    elseif (strcmpi(Type, 'full'))
        tree = KDTreeSearcher(GM.V');
        pointOnSourceIndex = tree.knnsearch(selectedPoint');
        pointOnSource = GM.V(:,pointOnSourceIndex);
    end
    clear GM;
    h = plot3(findobj('Tag','Target'), pointOnSource(1,:),...
        pointOnSource(2,:), pointOnSource(3,:), 'r.', 'MarkerSize', MarkerSize);
    set(h,'Tag','pt');
    if ~isempty(findobj('Tag','PropagTarget'))
        h = plot3(findobj('Tag','PropagTarget'), pointOnSource(1,:),...
            pointOnSource(2,:), pointOnSource(3,:), 'g.', 'MarkerSize', MarkerSize);
        set(h,'Tag','pt');
    end
    
    GM = mesh_list{1};
    if (strcmpi(Type, 'sample'))
        tree = KDTreeSearcher(GM.V(:,GM.Aux.VertSampInd)');
        pointOnTargetIndex = map(pointOnSourceIndex);
        pointOnTarget = GM.V(:,pointOnTargetIndex);
        pointOnTargetSampleIndex = tree.knnsearch(pointOnTarget');
        pointOnTargetIndex = GM.Aux.VertSampInd(pointOnTargetSampleIndex);
        pointOnTarget = GM.V(:,pointOnTargetIndex);
    elseif (strcmpi(Type, 'full'))
        pointOnTarget = GM.V(:,map(pointOnSourceIndex));
    end
    clear GM;
    h = plot3(findobj('Tag','Source'), pointOnTarget(1,:),...
        pointOnTarget(2,:), pointOnTarget(3,:), 'g.', 'MarkerSize', MarkerSize);
    set(h,'Tag','pt');
    if ~isempty(findobj('Tag','PropagSource'))
        h = plot3(findobj('Tag','PropagSource'), pointOnTarget(1,:),...
            pointOnTarget(2,:), pointOnTarget(3,:), 'g.', 'MarkerSize', MarkerSize);
        set(h,'Tag','pt');
    end
    
end

end



