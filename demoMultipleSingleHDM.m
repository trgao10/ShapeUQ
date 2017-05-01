clearvars;
% close all;
path(pathdef);
path(path, genpath('./utils'));

MapType = 'cPMST';
FeatureFix = '';
GroupLevel = 'Genus';
numLmk = 20;
dispLayout = [5,10];
GroupNames = {'Alouatta','Ateles','Brachyteles','Callicebus','Saimiri'};
% GroupNames = {'Alouatta'};


base_path = [pwd '/'];
data_path = '../DATA/HDM/';
sample_path = [data_path 'samples/'];
meshes_path = [data_path 'meshes/'];
spreadsheet_path = [data_path 'ClassificationTable.xlsx'];

taxa_file = [data_path 'hdm_taxa_table.mat'];
taxa_code = load(taxa_file);
taxa_code = taxa_code.taxa_code;

[~,ClTable,~] = xlsread(spreadsheet_path);
Names = {};
NamesByGroup = cell(1,length(GroupNames));
TaxaByGroup = cell(1,length(GroupNames));
for j=1:length(GroupNames)
    NamesJ = ClTable(strcmpi(ClTable(1:end,strcmpi(ClTable(1,:),GroupLevel)),GroupNames{j}),1);
    Names = [Names,NamesJ{:}];
    NamesByGroup{j} = NamesJ;
    TaxaByGroup{j} = cellfun(@(x) find(strcmpi(taxa_code, x)), NamesJ);
end

GroupSize = length(Names);
DiffMatrixSizeList = zeros(GroupSize,1);
TAXAinds = zeros(GroupSize,1);
NamesDelimit = zeros(GroupSize+1,2);
meshList = cell(1,GroupSize);
for j=1:GroupSize
    TAXAinds(j) = find(strcmpi(taxa_code,Names{j}));
    load([sample_path taxa_code{strcmpi(taxa_code,Names{j})} '.mat']);
    meshList{j} = G;
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

GroupSize = length(Names);

TAXAinds = zeros(GroupSize,1);
for j=1:GroupSize
    TAXAinds(j) = find(strcmpi(taxa_code,Names{j}));
end
MeshNames = taxa_code(TAXAinds); % match upper/lower cases

rigid_motions = load([data_path 'cPMSTinitRms.mat']);
R = rigid_motions.R(TAXAinds,TAXAinds);

figure;
h = zeros(1,prod(dispLayout));
for j=1:length(MeshNames)
    MeshName = MeshNames{j};
    G = Mesh('off',[meshes_path MeshName '.off']);    
    [GPLmkIdx,ptuq] = G.GetGPLmk(numLmk);
    
    h(j) = subplot(dispLayout(1),dispLayout(2),j);
    G.V = R{1,j} * G.V;
    G.ViewFunctionOnMesh(ptuq,struct('mode','native'));
    hold on
    scatter3(G.V(1,GPLmkIdx),G.V(2,GPLmkIdx),G.V(3,GPLmkIdx),20,'g','filled');
end

Link = linkprop(h, {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'CameraViewAngle'});
setappdata(gcf, 'StoreTheLink', Link);
