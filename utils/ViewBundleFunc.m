function [] = ViewBundleFunc(Names,ConsistScore,options)

sample_path = options.sample_path;
DisplayLayout = options.DisplayLayout;
GroupSize = length(Names);
mesh_list = cell(size(Names));
R = getoptions(options, 'R', repmat({eye(3)}, 1, GroupSize));
linkCamera = getoptions(options,'linkCamera','on');
DisplayOrient = getoptions(options,'DisplayOrient','Horizontal');

switch DisplayOrient
    case 'Vertical'
        DisplayOrder = reshape(1:DisplayLayout(1)*DisplayLayout(2), DisplayLayout(2), DisplayLayout(1));
        DisplayOrder = DisplayOrder';
        DisplayOrder = DisplayOrder(:);
    case 'Horizontal'
        DisplayOrder = 1:DisplayLayout(1)*DisplayLayout(2);
end

for i=1:GroupSize
    GM = load([sample_path Names{i} '.mat']);
    GM = GM.G;
    GM.V = R{i}*GM.V;
    mesh_list{i} = GM;
%     % align every tooth to the first one on the list
%     if (i==1)
%         mesh_list{i} = GM;
%     else
%         GM.V = R{i}*GM.V;
% %         GM.V = R{1,i}*GM.V;
%         mesh_list{i} = GM;
%     end
end

if (~isempty(findobj('Tag','BundleFunc')))
    camUpVector = get(gca, 'CameraUpVector');
    camPosition = get(gca, 'CameraPosition');
    camTarget = get(gca, 'CameraTarget');
    camViewAngle = get(gca, 'CameraViewAngle');
    figurePosition = get(gcf, 'Position');
else
    figurePosition = [10, 10, 800, 800];
end

figure('Unit', 'pixel', 'Position', figurePosition);
set(gcf, 'ToolBar', 'none');
h = zeros(size(mesh_list));

BlockShift = zeros(1,GroupSize);
for j=1:GroupSize
    BlockShift(j) = mesh_list{j}.nV;
end
BlockShift = cumsum(BlockShift);
BlockShift = [0 BlockShift];

for i=1:GroupSize
    color_data = ConsistScore((BlockShift(i)+1):BlockShift(i+1));
    h(i) = subplot(DisplayLayout(1), DisplayLayout(2), DisplayOrder(i));
    mesh_list{i}.draw(struct('FaceColor', 'interp', 'FaceVertexCData', color_data, 'CDataMapping','scaled', 'EdgeColor', 'none', 'FaceAlpha',1,'AmbientStrength',0.3,'SpecularStrength',0.0));
    hold on;
%     colormap jet(256);
    camlight('headlight');
    camlight(180,0);
    if strcmpi(options.names,'on')
        title(mesh_list{i}.Aux.name);
    end
end

if strcmpi(linkCamera, 'on')
    Link = linkprop(h, {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'CameraViewAngle'});
    setappdata(gcf, 'StoreTheLink', Link);
end

if (exist('camUpVector', 'var'))
    set(gca, 'CameraUpVector', camUpVector);
    set(gca, 'CameraPosition', camPosition);
    set(gca, 'CameraTarget', camTarget);
    set(gca, 'CameraViewAngle', camViewAngle);
    set(gcf, 'Tag', 'BundleFunc');
end

end


