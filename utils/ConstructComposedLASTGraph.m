function [ST] = ConstructComposedLASTGraph(DistMatrix,MapMatrix,TaxaCode,options)
%CONSTRUCTCPLASTGRAPH Summary of this function goes here
%   Detailed explanation goes here

DistMatrix = sparse(DistMatrix);
TrilDistMatrix = tril(DistMatrix, -1);
alpha = getoptions(options,'alpha',8);

OPT_WEIGHT = Inf;
RowNodes = 1:size(DistMatrix,1);
Weights = zeros(size(RowNodes));
for RootNode=1:size(DistMatrix,1)
    progressbar(RootNode,size(DistMatrix,1),20);
    
    GM = load([options.SamplePath TaxaCode{RootNode} '.mat']);
    GM = GM.G;
    ST = graphminspantree(TrilDistMatrix, RootNode, 'Method', 'Kruskal');
    [~,MinPaths,~] = graphshortestpath(ST,RootNode,'directed',false);
    DistRootJ = zeros(1,size(DistMatrix,1));
    for j=1:length(DistRootJ)
        GN = load([options.SamplePath TaxaCode{j} '.mat']);
        GN = GN.G;
        DistRootJ(j) = MapToDist(GM.V,GN.V,ComposeMapsAlongPath(MinPaths{j},MapMatrix),GM.Aux.VertArea);
    end
    if ~isnumeric(alpha)
        PointwiseDistortion = DistRootJ./DistMatrix(RootNode,:);
        alpha = mean(PointwiseDistortion);
    end
    HighDistortionInds = find(DistRootJ>DistMatrix(RootNode,:)*alpha);
    ST(RootNode, HighDistortionInds) = DistMatrix(RootNode, HighDistortionInds);
    ST = max(ST,ST');
    [~,~,PRED] = graphshortestpath(ST,RootNode,'directed',false);
    for j=1:length(RowNodes)
        if PRED(j)==0
            Weights(j) = 0;
        else
            Weights(j) = DistMatrix(RowNodes(j),PRED(j));
        end
    end
    if sum(Weights)<OPT_WEIGHT
        OPT_WEIGHT = sum(Weights);
        OPT_RootNode = RootNode;
        OPT_PRED = PRED;
        OPT_Weights = Weights;
    end
end
RowNodes(OPT_RootNode) = [];
OPT_PRED(OPT_RootNode) = [];
OPT_Weights(OPT_RootNode) = [];
SPT = sparse(RowNodes,OPT_PRED,OPT_Weights,size(DistMatrix,1),size(DistMatrix,2));
ST = max(SPT,SPT');
ST = tril(ST, -1);

end

