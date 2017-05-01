%% preparation
close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

%% setup parameters
BaseEps = 0.03;
BNN = 5;
FibrEps = 1e-3;
MapType = 'cPMST';
FeatureFix = 'Off';
numLmk = 7;
GroupLevel = 'Genus';
% GroupNames = {'Purgatorius','Pronothodectes','Tupaia','Lemur',...
%     'Microcebus','Cantius','Arctocebus','Adapis','Lepilemur',...
%     'Eosimias','Cynocephalus','Leptacodon','Nycticebus'};
% GroupNames = {'Euprimates','Primates','Dermoptera','Scandentia','Incertae sedis'};
% GroupNames = {'Purgatorius'};
% GroupNames = {'Purgatorius','Pronothodectes'};
GroupNames = {'Purgatorius','Pronothodectes','Tupaia','Lemur'};
% GroupNames = {'Purgatorius','Pronothodectes','Tupaia','Lemur',...
%     'Microcebus','Cantius','Arctocebus','Adapis','Lepilemur',...
%     'Eosimias','Cynocephalus'};
% GroupNames = {'Donrussellia','Cheirogaleus','Avahi','Eulemur',...
%     'Hapalemur','Loris','Nycticebus','Leptacodon'};
% GroupNames = {'Tupaia','Galago'};
% GroupNames = {'Purgatorius','Tupaia','Pronothodectes','Varecia','Microcebus','Lemur'};

%% setup paths
base_path = [pwd '/'];
data_path = '../DATA/PNAS/';
spreadsheet_path = [data_path 'ClassificationTable.xlsx'];
% sample_path = '../cPdist/samples/PNAS/';
sample_path = [data_path 'samples/'];
result_path = ['/media/trgao10/Storage/Work/MATLAB/ArchivedResults/PNAS/' MapType '/' 'FeatureFix' FeatureFix '/'];
soften_path = [result_path 'soften/'];

%% options that control the diffusion eigenvector visualization
options.sample_path = sample_path;
options.DisplayLayout = [4,4];
options.DisplayOrient = 'Horizontal';
options.boundary = 'on';
options.names = 'off';

%% load taxa codes
taxa_file = [data_path 'teeth_taxa_table.mat'];
taxa_code = load(taxa_file);
taxa_code = taxa_code.taxa_code;
GroupSize = length(taxa_code);
ChunkSize = 55; %% PNAS

%% useful inline functions
flat = @(x) x(:);
ChunkIdx = @(TAXAind1,TAXAind2) ceil(((TAXAind1-1)*GroupSize+TAXAind2)/ChunkSize);

%% parse GroupNames
[~,ClTable,~] = xlsread(spreadsheet_path);
Names = {};
NamesByGroup = cell(1,length(GroupNames));
for j=1:length(GroupNames)
    NamesJ = ClTable(strcmpi(ClTable(1:end,strcmpi(ClTable(1,:),GroupLevel)),GroupNames{j}),1);
    Names = [Names,NamesJ{:}];
    NamesByGroup{j} = NamesJ;
end

GroupSize = length(Names);
DiffMatrixSizeList = zeros(GroupSize,1);
TAXAinds = zeros(GroupSize,1);
NamesDelimit = zeros(GroupSize+1,2);
for j=1:GroupSize
    TAXAinds(j) = find(strcmpi(taxa_code,Names{j}));
    load([sample_path taxa_code{strcmpi(taxa_code,Names{j})} '.mat']);
    DiffMatrixSizeList(j) = G.nV;
    NamesDelimit(j+1,1) = NamesDelimit(j,2)+1;
    NamesDelimit(j+1,2) = NamesDelimit(j+1,1)+G.nV-1;
end
Names = taxa_code(TAXAinds); % match upper/lower cases
NamesDelimit(1,:) = [];
nVList = DiffMatrixSizeList;
nVListCumsum = cumsum(nVList);

PerGroupSize = zeros(1,length(GroupNames));
for j=1:length(NamesByGroup)
    for k=1:length(NamesByGroup{j})
        NamesByGroup{j}{k} = taxa_code{strcmpi(taxa_code,NamesByGroup{j}{k})};
    end
    PerGroupSize(j) = length(NamesByGroup{j});
end
CumsumPerGroupSize = cumsum(PerGroupSize);

%% collection rigid motions
rigid_motions = load([data_path 'cPMSTinitRms.mat']);
options.R = rigid_motions.R(TAXAinds(1),TAXAinds);

%% process base diffusion
load([result_path MapType 'DistMatrix.mat']);
if strcmpi(MapType,'cP')
    eval(['BaseDistMatrix = ' MapType 'DistMatrix(TAXAinds,TAXAinds);']);
else
    eval(['BaseDistMatrix = ImprDistMatrix(TAXAinds,TAXAinds);']);
end
BaseDistMatrix = BaseDistMatrix-diag(diag(BaseDistMatrix));

%%% only connect BNN-nearest-neighbors
[sDists,rowNNs] = sort(BaseDistMatrix,2);
sDists = sDists(:,2:(1+BNN));
rowNNs = rowNNs(:,2:(1+BNN));
BaseWeights = sparse(repmat((1:GroupSize)',1,BNN),rowNNs,sDists,GroupSize,GroupSize);
BaseWeights = min(BaseWeights, BaseWeights');
for j=1:GroupSize
    sDists(j,:) = BaseWeights(j,rowNNs(j,:));
end
sDists = exp(-sDists.^2/BaseEps);

%% build diffusion kernel matrix
DiffMatrixSize = sum(DiffMatrixSizeList);
DiffMatrixSizeList = cumsum(DiffMatrixSizeList);
DiffMatrixSizeList = [0; DiffMatrixSizeList];
GroupDelimit = zeros(length(GroupNames)+1,2);
for j=2:(length(GroupNames)+1)
    GroupDelimit(j,1) = GroupDelimit(j-1,2)+1;
    GroupDelimit(j,2) = DiffMatrixSizeList(CumsumPerGroupSize(j-1)+1);
end
GroupDelimit(1,:) = [];
DiffMatrixSizeList(end) = []; % treated as block shifts
DiffMatrixRowIdx = [];
DiffMatrixColIdx = [];
DiffMatrixVal = [];

VertAreaMeas = zeros(DiffMatrixSize,1);
meshList = cell(1,GroupSize);

cback = 0;
for j=1:GroupSize
    G1 = load([sample_path taxa_code{strcmpi(taxa_code,Names{j})} '.mat']);
    G1 = G1.G;
    meshList{j} = G1;
    meshList{j}.V = options.R{j}*meshList{j}.V;
    [Cgauss,Cmean] = G1.ComputeCurvature();
    VertAreaMeas(NamesDelimit(j,1):NamesDelimit(j,2)) = G1.Aux.VertArea'+abs(Cgauss)/sum(abs(Cgauss))+abs(Cmean)/sum(abs(Cmean));
    
    for nns = 1:BNN
        if (sDists(j,nns) == 0)
            continue;
        end
        k = rowNNs(j,nns);
        G2 = load([sample_path taxa_code{strcmpi(taxa_code,Names{k})} '.mat']); G2 = G2.G;
        
        %%% load texture coordinates
        TAXAind1 = TAXAinds(j);
        TAXAind2 = TAXAinds(k);
        load([soften_path 'soften_mat_' num2str(ChunkIdx(TAXAind1, TAXAind2)) '.mat']);
        AugKernel12 = cPSoftMapsMatrix{TAXAind1, TAXAind2};
        
        [rowIdx, colIdx, val] = find(AugKernel12);
        DiffMatrixRowIdx = [DiffMatrixRowIdx; rowIdx+DiffMatrixSizeList(j)];
        DiffMatrixColIdx = [DiffMatrixColIdx; colIdx+DiffMatrixSizeList(k)];
        DiffMatrixVal = [DiffMatrixVal; sDists(j,nns)*val];
        [rowIdx, colIdx, val] = find(AugKernel12');
        DiffMatrixRowIdx = [DiffMatrixRowIdx; rowIdx+DiffMatrixSizeList(k)];
        DiffMatrixColIdx = [DiffMatrixColIdx; colIdx+DiffMatrixSizeList(j)];
        DiffMatrixVal = [DiffMatrixVal; sDists(j,nns)*val];
    end
    
    for cc=1:cback
        fprintf('\b');
    end
    cback = fprintf(['%4d/' num2str(GroupSize) ' done.\n'],j);
end

H = sparse(DiffMatrixRowIdx,DiffMatrixColIdx,DiffMatrixVal,DiffMatrixSize,DiffMatrixSize);
clear DiffMatrixColIdx DiffMatrixRowIdx DiffMatrixVal rowIdx colIdx val
clear TextureCoords1Matrix TextureCoords2Matrix

%% eigen-decomposition
fullPhi = H;
sqrtD = sparse(1:DiffMatrixSize,1:DiffMatrixSize,sqrt(sum(H)));
invD = sparse(1:DiffMatrixSize,1:DiffMatrixSize,1./sum(H));
sqrtInvD = sparse(1:DiffMatrixSize,1:DiffMatrixSize,1./sqrt(sum(H)));
H = sqrtInvD*H*sqrtInvD;
H = (H+H')/2;

eigopt = struct('isreal',1,'issym',1,'maxit',5000,'disp',0);
tic;
[U, lambda] = eigs(H, 101, 'LM', eigopt);
lambda = diag(lambda);
disp(['Eigs completed in ' num2str(toc) ' seconds']);

%% construct HDM
sqrtInvD(isinf(sqrtInvD)) = 0;
BundleHDM = sqrtInvD*U(:,2:15);
kdtreeList = cell(1,GroupSize);
for j=1:GroupSize
    kdtreeList{j} = kdtree_build(BundleHDM(NamesDelimit(j,1):NamesDelimit(j,2),:));
end

%% calculate pointwise uncertainty over all shapes in the collection
VertAreaMeas = sparse(1:DiffMatrixSize,1:DiffMatrixSize,VertAreaMeas,DiffMatrixSize,DiffMatrixSize);
ptuq = sum(fullPhi.*(VertAreaMeas*fullPhi))';
% ptuq = sum(H.*(sparse(1:DiffMatrixSize,1:DiffMatrixSize,VertAreaMeas,DiffMatrixSize,DiffMatrixSize)*H))';
%%%% locate on which shape the largest ptuq occurs, pick as InitLmkIdx
[~,InitLmkIdx] = max(ptuq);
InitLmkShapeIdx = find((NamesDelimit(:,1)<=InitLmkIdx) & (NamesDelimit(:,2)>=InitLmkIdx));
AllInitLmkIdx = zeros(GroupSize,1);
AllInitLmkIdx(InitLmkShapeIdx) = InitLmkIdx-NamesDelimit(InitLmkShapeIdx,1)+1;
for j=1:GroupSize
    if (j == InitLmkShapeIdx)
        continue;
    end
    AllInitLmkIdx(j) = kdtree_nearest_neighbor(kdtreeList{j}, BundleHDM(InitLmkIdx,:));
end

%%%% plot the peak location
% h = zeros(1,GroupSize);
% figure;
% for j=1:GroupSize
%     h(j) = subplot(options.DisplayLayout(1),options.DisplayLayout(2),j);
%     meshList{j}.ViewFunctionOnMesh(full(ptuq(NamesDelimit(j,1):NamesDelimit(j,2))),struct('mode','native'));
%     hold on
%     localInitLmkIdx = AllInitLmkIdx(j);
%     scatter3(meshList{j}.V(1,localInitLmkIdx),...
%         meshList{j}.V(2,localInitLmkIdx),...
%         meshList{j}.V(3,localInitLmkIdx),...
%         20,'g','filled');
% end
% Link = linkprop(h, {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'CameraViewAngle'});
% setappdata(gcf, 'StoreTheLink', Link);

GPLmkIdx = zeros(GroupSize,numLmk);
GPLmkIdx(:,1) = AllInitLmkIdx+NamesDelimit(:,1)-1;
locGPLmkIdx = zeros(GroupSize,numLmk);
locGPLmkIdx(:,1) = AllInitLmkIdx;

cback = 0;
for k=2:numLmk
    for cc=1:cback
        fprintf('\b');
    end
    cback = fprintf('Landmark: %4d/%d\n',k,numLmk);
    
    lmkPhi = fullPhi(:,flat(GPLmkIdx(:,1:(k-1))));
    projNullPhi = VertAreaMeas-VertAreaMeas*lmkPhi*((lmkPhi'*VertAreaMeas*lmkPhi)\(lmkPhi'*VertAreaMeas));
    ptuq = sum(fullPhi.*(projNullPhi*fullPhi))';
    
    [~,MaxLmkIdx] = max(ptuq);
    MaxLmkShapeIdx = find((NamesDelimit(:,1)<=MaxLmkIdx) & (NamesDelimit(:,2)>=MaxLmkIdx));
    AllInitLmkIdx = zeros(GroupSize,1);
    AllInitLmkIdx(MaxLmkShapeIdx) = MaxLmkIdx-NamesDelimit(MaxLmkShapeIdx,1)+1;
    for j=1:GroupSize
        if (j == MaxLmkShapeIdx)
            continue;
        end
        AllInitLmkIdx(j) = kdtree_nearest_neighbor(kdtreeList{j}, BundleHDM(MaxLmkIdx,:));
    end
    
    GPLmkIdx(:,k) = AllInitLmkIdx+NamesDelimit(:,1)-1;
    locGPLmkIdx(:,k) = AllInitLmkIdx;
    
%     h = zeros(1,GroupSize);
%     figure('Name',['Landmark ' num2str(k)]);
%     for j=1:GroupSize
%         h(j) = subplot(options.DisplayLayout(1),options.DisplayLayout(2),j);
%         meshList{j}.ViewFunctionOnMesh(full(ptuq(NamesDelimit(j,1):NamesDelimit(j,2))),struct('mode','native'));
%         hold on
%         localInitLmkIdx = AllInitLmkIdx(j);
%         scatter3(meshList{j}.V(1,localInitLmkIdx),...
%             meshList{j}.V(2,localInitLmkIdx),...
%             meshList{j}.V(3,localInitLmkIdx),...
%             20,'g','filled');
%     end
%     Link = linkprop(h, {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'CameraViewAngle'});
%     setappdata(gcf, 'StoreTheLink', Link);
end

%% plot extracted landmarks
h = zeros(1,GroupSize);
figure('Name',['Landmark ' num2str(k)]);
colorsList = [1,0,0;0,1,0;0,0,1;1,1,0;1,0,1;0,1,1;1,1,1];
colorsList = [colorsList;colorsList/2];
for j=1:GroupSize
    h(j) = subplot(options.DisplayLayout(1),options.DisplayLayout(2),j);
    color_data = repmat([0.9, 0.9, 0.8], meshList{j}.nV, 1);
    meshList{j}.draw(struct('FaceColor', 'interp',...
        'FaceVertexCData', color_data, 'CDataMapping','scaled',...
        'EdgeColor', 'none', 'FaceAlpha', 1,...
        'AmbientStrength',0.3,'SpecularStrength',0.0));
    lighting phong
    camlight('headlight');
    camlight(180,0);
    
    hold on
    for k=1:numLmk
        scatter3(meshList{j}.V(1,locGPLmkIdx(j,k)),...
            meshList{j}.V(2,locGPLmkIdx(j,k)),...
            meshList{j}.V(3,locGPLmkIdx(j,k)),...
            30,colorsList(k,:),'filled');
    end
end
Link = linkprop(h, {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'CameraViewAngle'});
setappdata(gcf, 'StoreTheLink', Link);


