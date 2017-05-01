clearvars;
% close all;
path(pathdef);
path(path, genpath('./utils'));

data_path = '~/Work/MATLAB/DATA/PNAS/';
sample_path = [data_path 'samples/'];
meshes_path = [data_path 'meshes/'];
% MeshName = 'D09';
% MeshName = 'B03';
% MeshName = 'x23';
MeshName = 'h08';

load([sample_path MeshName '.mat']);
numLmk = 50;
% lambda = 0.5;

figure;

%% extract GP landmarks
% [GPLmkIdx,ptuq] = G.GetGPLmk(numLmk,lambda);
[GPLmkIdx,ptuq] = G.GetGPLmk(numLmk);

h(1) = subplot(1,2,1);
G.ViewFunctionOnMesh(full(ptuq),struct('mode','native'));
hold on
scatter3(G.V(1,GPLmkIdx(1)),G.V(2,GPLmkIdx(1)),G.V(3,GPLmkIdx(1)),20,'r','filled');
scatter3(G.V(1,GPLmkIdx(2:end)),G.V(2,GPLmkIdx(2:end)),G.V(3,GPLmkIdx(2:end)),20,'g','filled');
title('Gaussian Process', 'Interpreter', 'Latex', 'FontSize', 15);

%% compare with FPS landmarks
GFPSIdx = G.GeodesicFarthestPointSampling(numLmk,GPLmkIdx(1));

h(2) = subplot(1,2,2);
G.ViewFunctionOnMesh(full(ptuq),struct('mode','native'));
hold on
scatter3(G.V(1,GFPSIdx(1)),G.V(2,GFPSIdx(1)),G.V(3,GFPSIdx(1)),20,'r','filled');
scatter3(G.V(1,GFPSIdx(2:end)),G.V(2,GFPSIdx(2:end)),G.V(3,GFPSIdx(2:end)),20,'g','filled');
title('Geodesic Farthese Point Sampling', 'Interpreter', 'Latex', 'FontSize', 15);

%% sync the two subplots
Link = linkprop(h, {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'CameraViewAngle'});
setappdata(gcf, 'StoreTheLink', Link);
