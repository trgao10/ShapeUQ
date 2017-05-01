%% preparation
clear vars;
% close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

data_path = '../data/PNAS/';
sample_path = [data_path 'samples/'];
meshes_path = [data_path 'meshes/'];

%% set parameters
numGPLmk = 50;
WKSNN = 5; %%% number of closest WKS descriptors for each landmark
Names = {'w09', 'w10'};
options = struct('FeatureType','ConfMax',...
    'NumDensityPnts', 100,...
    'AngleIncrement', 0.05,...
    'NumFeatureMatch', 4,...
    'GaussMinMatch', 'off');

%% load meshes (already flattend --- see field "Aux.UniformizationV")
Gs = cell(2,1);
for j=1:2
    load([sample_path Names{j} '.mat']);
    Gs{j} = G;
end

%% show unaligned uniformization results
figure;
for j=1:2
    subplot(1,2,j);
    D = Mesh('VF',Gs{j}.Aux.UniformizationV,Gs{j}.F);
    D.draw();
    hold on
    scatter3(D.V(1,Gs{j}.Aux.ConfMaxInds),...
        D.V(2,Gs{j}.Aux.ConfMaxInds),...
        D.V(3,Gs{j}.Aux.ConfMaxInds),...
        20,'g','filled');
end

%% get initial alignments
rslt12 = Gs{1}.ComputeContinuousProcrustes(Gs{2},options);
TextureCoords = {rslt12.TextureCoords1,rslt12.TextureCoords2};

rsltMeshes = cell(1,2);
figure;
for j=1:2
    subplot(1,2,j);
    rsltMeshes{j} = Mesh('VF',[TextureCoords{j};zeros(1,Gs{j}.nV)],Gs{j}.F);
    rsltMeshes{j}.draw();
    hold on
    scatter3(rsltMeshes{j}.V(1,Gs{j}.Aux.ConfMaxInds),...
        rsltMeshes{j}.V(2,Gs{j}.Aux.ConfMaxInds),...
        rsltMeshes{j}.V(3,Gs{j}.Aux.ConfMaxInds),...
        20,'g','filled');
end

%% try deforming the boundaries into perfect circles using CPMS
%%% can't use this uniformization as input to CPD computation because it
%%% will ruin the maps output from CPD for stupid reasons

GoodBoundaryMeshes = cell(1,2);
figure;
for j=1:2
    subplot(1,2,j);
    GoodBoundaryMeshes{j} = Mesh('VF',[TextureCoords{j};zeros(1,Gs{j}.nV)],Gs{j}.F).ComputeCPMS().ComputeUniformization(struct('method','LSCM','boundary_conditions','disc'));
    FlatVertices = GoodBoundaryMeshes{j}.V;
    [U,~,V] = svd(FlatVertices*rsltMeshes{j}.V');
    RJ = V*U';
    FlatVertices = RJ*FlatVertices;
    GoodBoundaryMeshes{j}.V = FlatVertices;
    
    GoodBoundaryMeshes{j}.draw();
    hold on
    scatter3(GoodBoundaryMeshes{j}.V(1,Gs{j}.Aux.ConfMaxInds),...
        GoodBoundaryMeshes{j}.V(2,Gs{j}.Aux.ConfMaxInds),...
        GoodBoundaryMeshes{j}.V(3,Gs{j}.Aux.ConfMaxInds),...
        20,'g','filled');
end

%% extract an equal number of GP landmarks on each mesh, together with WKS
for j=1:2
    Gs{j}.Aux.GPLmkIdx = Gs{j}.GetGPLmk(numLmk);
    Gs{j}.Aux.WKS = Gs{j}.ComputeWKS([]);
end

%% find WSKNN nearest neighors for each landmark on Gs{1}
distMat = pdist2(Gs{1}.Aux.WKS(Gs{1}.Aux.GPLmkIdx,:),...
                 Gs{2}.Aux.WKS(Gs{2}.Aux.GPLmkIdx,:));
[sDists,rowNNs] = sort(distMat,2);
%%% the $k$-th row of corrMap are the indices of the WKSNN-closest
%%% candidate landmarks on Gs{2} that will be potentially matches with the
%%% $k$-th landmark on Gs{1}
corrMap = rowNNs(:,1:WKSNN);

%% randomly pick a landmark on Gs{1} to see if the WKS matching makes sense
%%% (At least one of the green points on the right panel should match
%%% reasonably well with the red points on the left panel; keep in mind
%%% that the two meshes are not aligned with each other so a mental "flip"
%%% is probably needed.)
randLmkIdx = randi(numGPLmk,1);

figure;
h = zeros(2,1);
h(1) = subplot(1,2,1);
Gs{1}.draw();
hold on
scatter3(Gs{1}.V(1,Gs{1}.Aux.GPLmkIdx(randLmkIdx)),...
    Gs{1}.V(2,Gs{1}.Aux.GPLmkIdx(randLmkIdx)),...
    Gs{1}.V(3,Gs{1}.Aux.GPLmkIdx(randLmkIdx)),...
    20,'r','filled');

h(2) = subplot(1,2,2);
Gs{2}.draw();
hold on
scatter3(Gs{2}.V(1,Gs{2}.Aux.GPLmkIdx(corrMap(randLmkIdx,:))),...
    Gs{2}.V(2,Gs{2}.Aux.GPLmkIdx(corrMap(randLmkIdx,:))),...
    Gs{2}.V(3,Gs{2}.Aux.GPLmkIdx(corrMap(randLmkIdx,:))),...
    20,'g','filled');

%% apply l0-feature matching with bounded distortion maps
%%% refer to Gs{1 | 2}.Aux.GPLmkIdx for candidate feature points
%%% read off candidate matchings from the matrix "corrMap"
%%% work with GoodBoundaryMeshes to perform in-plane BD matching


