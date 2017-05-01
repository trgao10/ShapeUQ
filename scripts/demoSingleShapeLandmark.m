clearvars;
close all;
path(pathdef);
path(path, genpath('./utils'));

% flat = @(x) x(:);

data_path = '~/Work/MATLAB/DATA/PNAS/';
sample_path = [data_path 'samples/'];
meshes_path = [data_path 'meshes/'];
% samples_path = '~/Work/MATLAB/DATA/PNAS/samples/';
% MeshName = 'B03';
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% use D09 to make the point that these landmakrs are not simply "pooling"
%%% together all manually crafted geometric features
MeshName = 'D09'; 
% MeshName = 'a19';
load([sample_path MeshName '.mat']);
numLmk = 18;
G.Aux.LMK = G.GetLandmarksFromPNAS([data_path 'landmarks_teeth.mat'],[meshes_path G.Aux.name '_sas.off'],numLmk);

[Cgauss,Cmean,Umin,Umax,Cmin,Cmax,Normal] = G.ComputeCurvature();
DNE = Cmin.^2+Cmax.^2;
DNE(DNE<median(DNE)) = 0;

%%% check uncertainty of another landmark
% [~,TriArea] = G.ComputeSurfaceArea();
% VertArea = G.F2V'*TriArea;
% VertAreaMeas = diag(G.Aux.VertArea).*(diag(abs(Cgauss)/sum(abs(Cgauss))+abs(Cmean)/sum(abs(Cmean))));
VertAreaMeas = diag(G.Aux.VertArea)+diag(abs(Cgauss)/sum(abs(Cgauss))+abs(Cmean)/sum(abs(Cmean)));
% VertAreaMeas = diag(G.Aux.VertArea)+diag(abs(Cgauss)/sum(abs(Cgauss))+abs(Cmean)/sum(abs(Cmean))+sqrt(DNE)/sum(sqrt(DNE)));
% VertAreaMeas = diag(G.Aux.VertArea)+diag(abs(Cgauss)/sum(abs(Cgauss))+abs(Cmean)/sum(abs(Cmean)));
% VertAreaMeas = diag(G.Aux.VertArea)+diag(abs(Cgauss)+abs(Cmean)+abs(G.Aux.Conf));
% VertAreaMeas = diag(G.Aux.VertArea)+diag(abs(Cgauss)+abs(Cmean)+abs(G.Aux.Conf)/sum(abs(G.Aux.Conf)));
% VertAreaMeas = diag(G.Aux.VertArea)+diag(abs(Cgauss)/sum(abs(Cgauss))+abs(Cmean)/sum(abs(Cmean))+abs(G.Aux.Conf)/sum(abs(G.Aux.Conf)));
% VertAreaMeas = diag(G.Aux.VertArea)+diag(DNE/max(DNE));

[EdgeIdxI,EdgeIdxJ] = find(tril(G.E));
bandwidth = mean(sqrt(sum((G.V(:,EdgeIdxI)-G.V(:,EdgeIdxJ)).^2)))/3;
% bandwidth = 0.01;
PDistMat = squareform(pdist(G.V'));
fullPhi = exp(-PDistMat.^2/bandwidth);

%%% pick initial landmark
% InitLmkIdx = G.Aux.ConfMaxInds;
% % InitLmkIdx = G.Aux.ConfMaxInds;

ptuq = sum(fullPhi.*(VertAreaMeas*fullPhi))';
[~,InitLmkIdx] = max(ptuq);

GPLmkIdx = InitLmkIdx;
inferLmkIdx = [];
cback = 0;
for j=1:(numLmk-length(InitLmkIdx))
    for cc=1:cback
        fprintf('\b');
    end
    cback = fprintf('Landmark: %4d\n',j);
    
    lmkPhi = fullPhi(:,GPLmkIdx);
    projNullPhi = VertAreaMeas-VertAreaMeas*lmkPhi*((lmkPhi'*VertAreaMeas*lmkPhi)\(lmkPhi'*VertAreaMeas));
    ptuq = sum(fullPhi.*(projNullPhi*fullPhi))';
    [~,maxUQIdx] = max(ptuq);
    GPLmkIdx = [GPLmkIdx;maxUQIdx];
    inferLmkIdx = [inferLmkIdx;maxUQIdx];
%     clf
%     G.ViewFunctionOnMesh(ptuq,struct('mode','native'));
%     hold on
%     scatter3(G.V(1,InitLmkIdx),G.V(2,InitLmkIdx),G.V(3,InitLmkIdx),20,'r','filled');
%     scatter3(G.V(1,inferLmkIdx),G.V(2,inferLmkIdx),G.V(3,inferLmkIdx),20,'g','filled');
%     pause
end

lmkPhi = fullPhi(:,GPLmkIdx);
projNullPhi = VertAreaMeas-VertAreaMeas*lmkPhi*((lmkPhi'*VertAreaMeas*lmkPhi)\(lmkPhi'*VertAreaMeas));
ptuq = sum(fullPhi.*(projNullPhi*fullPhi))';

%%
% clf
figure;
color_data = repmat([0.9, 0.9, 0.8], G.nV, 1);
h = zeros(1,6);

h(1) = subplot(2,3,1);
G.draw(struct('FaceColor', 'interp',...
    'FaceVertexCData', color_data, 'CDataMapping','scaled',...
    'EdgeColor', 'none', 'FaceAlpha', 1,...
    'AmbientStrength',0.3,'SpecularStrength',0.0));
hold on
scatter3(G.V(1,G.Aux.LMK),G.V(2,G.Aux.LMK),G.V(3,G.Aux.LMK),20,'k','filled');
% material shiny
lighting phong
set(gca, 'CameraUpVector',  [-0.2666,-0.9468, 0.1804]);
set(gca, 'CameraPosition',  [ 5.4366,-1.4825, 0.2070]);
set(gca, 'CameraTarget',    [ 0.0059, 0.0039,-0.0166]);
set(gca, 'CameraViewAngle', 10.0659);
camlight('headlight');
camlight(180,0);
title('Observer', 'Interpreter', 'Latex');

h(2) = subplot(2,3,2);
G.draw(struct('FaceColor', 'interp',...
    'FaceVertexCData', color_data, 'CDataMapping','scaled',...
    'EdgeColor', 'none', 'FaceAlpha', 1,...
    'AmbientStrength',0.3,'SpecularStrength',0.0));
hold on
scatter3(G.V(1,G.Aux.GaussMaxInds),G.V(2,G.Aux.GaussMaxInds),G.V(3,G.Aux.GaussMaxInds),20,'k','filled');
% material shiny
lighting phong
set(gca, 'CameraUpVector',  [-0.2666,-0.9468, 0.1804]);
set(gca, 'CameraPosition',  [ 5.4366,-1.4825, 0.2070]);
set(gca, 'CameraTarget',    [ 0.0059, 0.0039,-0.0166]);
set(gca, 'CameraViewAngle', 10.0659);
camlight('headlight');
camlight(180,0);
title('GaussMax', 'Interpreter', 'Latex');

h(3) = subplot(2,3,3);
G.draw(struct('FaceColor', 'interp',...
    'FaceVertexCData', color_data, 'CDataMapping','scaled',...
    'EdgeColor', 'none', 'FaceAlpha', 1,...
    'AmbientStrength',0.3,'SpecularStrength',0.0));
hold on
scatter3(G.V(1,G.Aux.GaussMinInds),G.V(2,G.Aux.GaussMinInds),G.V(3,G.Aux.GaussMinInds),20,'k','filled');
% material shiny
lighting phong
set(gca, 'CameraUpVector',  [-0.2666,-0.9468, 0.1804]);
set(gca, 'CameraPosition',  [ 5.4366,-1.4825, 0.2070]);
set(gca, 'CameraTarget',    [ 0.0059, 0.0039,-0.0166]);
set(gca, 'CameraViewAngle', 10.0659);
camlight('headlight');
camlight(180,0);
title('GaussMin', 'Interpreter', 'Latex');

h(4) = subplot(2,3,4);
G.draw(struct('FaceColor', 'interp',...
    'FaceVertexCData', color_data, 'CDataMapping','scaled',...
    'EdgeColor', 'none', 'FaceAlpha', 1,...
    'AmbientStrength',0.3,'SpecularStrength',0.0));
hold on
scatter3(G.V(1,GPLmkIdx),G.V(2,GPLmkIdx),G.V(3,GPLmkIdx),20,'k','filled');
% material shiny
lighting phong
set(gca, 'CameraUpVector',  [-0.2666,-0.9468, 0.1804]);
set(gca, 'CameraPosition',  [ 5.4366,-1.4825, 0.2070]);
set(gca, 'CameraTarget',    [ 0.0059, 0.0039,-0.0166]);
set(gca, 'CameraViewAngle', 10.0659);
camlight('headlight');
camlight(180,0);
title('GPLmk', 'Interpreter', 'Latex');

h(5) = subplot(2,3,5);
G.draw(struct('FaceColor', 'interp',...
    'FaceVertexCData', color_data, 'CDataMapping','scaled',...
    'EdgeColor', 'none', 'FaceAlpha', 1,...
    'AmbientStrength',0.3,'SpecularStrength',0.0));
hold on
scatter3(G.V(1,G.Aux.ConfMaxInds),G.V(2,G.Aux.ConfMaxInds),G.V(3,G.Aux.ConfMaxInds),20,'k','filled');
% material shiny
lighting phong
set(gca, 'CameraUpVector',  [-0.2666,-0.9468, 0.1804]);
set(gca, 'CameraPosition',  [ 5.4366,-1.4825, 0.2070]);
set(gca, 'CameraTarget',    [ 0.0059, 0.0039,-0.0166]);
set(gca, 'CameraViewAngle', 10.0659);
camlight('headlight');
camlight(180,0);
title('ConfMax', 'Interpreter', 'Latex');

h(6) = subplot(2,3,6);
G.draw(struct('FaceColor', 'interp',...
    'FaceVertexCData', color_data, 'CDataMapping','scaled',...
    'EdgeColor', 'none', 'FaceAlpha', 1,...
    'AmbientStrength',0.3,'SpecularStrength',0.0));
hold on
scatter3(G.V(1,G.Aux.ADMaxInds),G.V(2,G.Aux.ADMaxInds),G.V(3,G.Aux.ADMaxInds),20,'k','filled');
% material shiny
lighting phong
set(gca, 'CameraUpVector',  [-0.2666,-0.9468, 0.1804]);
set(gca, 'CameraPosition',  [ 5.4366,-1.4825, 0.2070]);
set(gca, 'CameraTarget',    [ 0.0059, 0.0039,-0.0166]);
set(gca, 'CameraViewAngle', 10.0659);
camlight('headlight');
camlight(180,0);
title('ADMax', 'Interpreter', 'Latex');

Link = linkprop(h, {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'CameraViewAngle'});
setappdata(gcf, 'StoreTheLink', Link);


figure;
G.ViewFunctionOnMesh(ptuq,struct('mode','native'));
hold on
scatter3(G.V(1,InitLmkIdx),G.V(2,InitLmkIdx),G.V(3,InitLmkIdx),20,'r','filled');
scatter3(G.V(1,inferLmkIdx),G.V(2,inferLmkIdx),G.V(3,inferLmkIdx),20,'g','filled');

% figure;
% G.draw();
% hold on
% scatter3(G.V(1,lmkIdx),G.V(2,lmkIdx),G.V(3,lmkIdx),20,'g','filled');
% G.ViewFunctionOnMesh(Cgauss,struct('mode','rb'));

