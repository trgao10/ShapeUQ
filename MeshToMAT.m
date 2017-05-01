clearvars;
% close all;
path(pathdef);
path(path, genpath('./utils'));

data_path = '~/Work/MATLAB/DATA/SCAPE/';
meshes_path = [data_path 'meshes/'];
samples_path = [data_path 'samples/'];

[MeshNames,suffix] = getFileNames(meshes_path);
for j=1:length(MeshNames)
    G = Mesh('off', [meshes_path MeshNames{j} suffix]);
    G.Centralize('ScaleArea');
    [~,TriArea] = G.ComputeSurfaceArea();
    G.Aux.VertArea = G.F2V'*TriArea;
    save([samples_path MeshNames{j} '.mat'], 'G');
end
