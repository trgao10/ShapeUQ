clearvars;
close all;
path(pathdef);
path(path, genpath('./utils'));

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
lambda = 0.5;
G.Aux.LMK = G.GetLandmarksFromPNAS([data_path 'landmarks_teeth.mat'],[meshes_path G.Aux.name '_sas.off'],numLmk);

[Cgauss,Cmean,Umin,Umax,Cmin,Cmax,Normal] = G.ComputeCurvature();

Lambda = diag(G.Aux.VertArea).*(lambda*diag(abs(Cgauss)/sum(abs(Cgauss))+(1-lambda)*abs(Cmean)/sum(abs(Cmean))));

[EdgeIdxI,EdgeIdxJ] = find(tril(G.E));
bandwidth = mean(sqrt(sum((G.V(:,EdgeIdxI)-G.V(:,EdgeIdxJ)).^2)))/2;
PDistMat = squareform(pdist(G.V'));
fullPhi = exp(-PDistMat.^2/bandwidth);

fullMatProd = fullPhi * Lambda * fullPhi;
KernelTrace = diag(fullMatProd);
% L = diag(sqrt(diag(Lambda))) * fullPhi;
GPLmkIdx = zeros(1,numLmk);
% lmkPhi = 0;

cback = 0;
for j=1:numLmk
    for cc=1:cback
        fprintf('\b');
    end
    cback = fprintf('Landmark: %4d\n',j);
    
%     lmkL = diag(sqrt(diag(Lambda))) * lmkPhi;
    if j == 1
        ptuq = KernelTrace;
    else
        ptuq = KernelTrace - sum(fullMatProd(GPLmkIdx(1:(j-1)),:)...
            .*(fullMatProd(GPLmkIdx(1:(j-1)),GPLmkIdx(1:(j-1)))\fullMatProd(GPLmkIdx(1:(j-1)),:)),1)';
    end
    [~,maxUQIdx] = max(ptuq);
    GPLmkIdx(j) = maxUQIdx;
    
%     lmkPhi = fullPhi(:,GPLmkIdx(1:j));
%     GPLmkIdx = [GPLmkIdx;maxUQIdx];
%     inferLmkIdx = [inferLmkIdx;maxUQIdx];
%     clf
%     G.ViewFunctionOnMesh(ptuq,struct('mode','native'));
%     hold on
%     scatter3(G.V(1,InitLmkIdx),G.V(2,InitLmkIdx),G.V(3,InitLmkIdx),20,'r','filled');
%     scatter3(G.V(1,inferLmkIdx),G.V(2,inferLmkIdx),G.V(3,inferLmkIdx),20,'g','filled');
%     pause
end

InitLmkIdx = GPLmkIdx(1);
inferLmkIdx = GPLmkIdx(2:end);

ptuq = KernelTrace - sum(fullMatProd(GPLmkIdx,:)...
    .*(fullMatProd(GPLmkIdx,GPLmkIdx)\fullMatProd(GPLmkIdx,:)),1)';

% lmkPhi = fullPhi(:,GPLmkIdx);
% projNullPhi = Lambda-Lambda*lmkPhi*((lmkPhi'*Lambda*lmkPhi)\(lmkPhi'*Lambda));
% ptuq = sum(fullPhi.*(projNullPhi*fullPhi))';

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

