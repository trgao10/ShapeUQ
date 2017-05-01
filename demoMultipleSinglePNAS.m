clearvars;
% close all;
path(pathdef);
path(path, genpath('./utils'));

data_path = '~/Work/MATLAB/DATA/PNAS/';
sample_path = [data_path 'samples/'];
meshes_path = [data_path 'meshes/'];
taxa_file = [data_path 'teeth_taxa_table.mat'];
taxa_code = load(taxa_file);
taxa_code = taxa_code.taxa_code;

numLmk = 16;
DisplayLayout = [2,4];
DisplayOrient = 'Horizontal';
% MeshNames = {'Q18','Q19','w01','w02'};
MeshNames = {'Q18','Q19','w01','w02',...
             'A15','A16','B19','B20'};

GroupSize = length(MeshNames);

TAXAinds = zeros(GroupSize,1);
for j=1:GroupSize
    TAXAinds(j) = find(strcmpi(taxa_code,MeshNames{j}));
end
MeshNames = taxa_code(TAXAinds); % match upper/lower cases

rigid_motions = load([data_path 'cPMSTinitRms.mat']);
RCell = rigid_motions.R;
for j=1:size(RCell,1)
    for k=1:size(RCell,2)
        RCell{k,j} = RCell{j,k}';
    end
end

[wGCL, wGCL_Dvec, wGCL_W] = assembleGCL(ones(size(RCell))-eye(size(RCell)), RCell, 3);
[RelaxSolCell_SPC, RelaxSolMat_SPC] = syncSpecRelax(wGCL, 3, wGCL_Dvec);

% R = rigid_motions.R(TAXAinds,TAXAinds);

switch DisplayOrient
    case 'Vertical'
        DisplayOrder = reshape(1:DisplayLayout(1)*DisplayLayout(2), DisplayLayout(2), DisplayLayout(1));
        DisplayOrder = DisplayOrder';
        DisplayOrder = DisplayOrder(:);
    case 'Horizontal'
        DisplayOrder = 1:DisplayLayout(1)*DisplayLayout(2);
end

figure;
h = zeros(1,prod(DisplayLayout));
for j=1:length(MeshNames)
    MeshName = MeshNames{j};
    G = Mesh('off',[meshes_path MeshName '_sas.off']);    
    [GPLmkIdx,ptuq] = G.GetGPLmk(numLmk);
    
    h(j) = subplot(DisplayLayout(1),DisplayLayout(2),DisplayOrder(j));
    G.V = RelaxSolCell_SPC{TAXAinds(j)} * G.V;
%     G.V = R{1,j} * G.V;
    G.ViewFunctionOnMesh(ptuq,struct('mode','native'));
    hold on
    scatter3(G.V(1,GPLmkIdx),G.V(2,GPLmkIdx),G.V(3,GPLmkIdx),20,'g','filled');
end

Link = linkprop(h, {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'CameraViewAngle'});
setappdata(gcf, 'StoreTheLink', Link);
