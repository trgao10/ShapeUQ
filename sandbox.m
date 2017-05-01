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
MeshNames = {'Q18','Q19','w01','w02',...
             'A15','A16','B19','B20'};

GroupSize = length(MeshNames);

TAXAinds = zeros(GroupSize,1);
for j=1:GroupSize
    TAXAinds(j) = find(strcmpi(taxa_code,MeshNames{j}));
end
MeshNames = taxa_code(TAXAinds); % match upper/lower cases

rigid_motions = load([data_path 'cPMSTinitRms.mat']);
R = rigid_motions.R(TAXAinds,TAXAinds);

figure;
dispLayout = [2,4];
h = zeros(1,prod(dispLayout));
for j=1:length(MeshNames)
    MeshName = MeshNames{j};
    G = Mesh('off',[meshes_path MeshName '_sas.off']);    
    [GPLmkIdx,ptuq] = G.GetGPLmk(numLmk);
    
    h(j) = subplot(dispLayout(1),dispLayout(2),j);
    G.V = R{1,j} * G.V;
    G.ViewFunctionOnMesh(ptuq,struct('mode','native'));
    hold on
    scatter3(G.V(1,GPLmkIdx),G.V(2,GPLmkIdx),G.V(3,GPLmkIdx),20,'g','filled');
end

Link = linkprop(h, {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'CameraViewAngle'});
setappdata(gcf, 'StoreTheLink', Link);


% data_path = '~/Work/MATLAB/DATA/AST_CALC/calc/';
% meshes_path = [data_path 'meshes/'];
% samples_path = [data_path 'samples/'];
% 
% [MeshNames,suffix] = getFileNames(meshes_path);
% for j=1:length(MeshNames)
%     G = Mesh('off', [meshes_path MeshNames{j} suffix]);
%     G.Centralize('ScaleArea');
%     [~,TriArea] = G.ComputeSurfaceArea();
%     G.Aux.VertArea = G.F2V'*TriArea;
%     save([samples_path MeshNames{j} '.mat'], 'G');
% end
