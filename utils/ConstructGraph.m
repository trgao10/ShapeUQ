function [ST,PRED] = ConstructGraph(DistMatrix,Type,options)
%CONSTRUCTGRAPH Summary of this function goes here
%   Detailed explanation goes here

if nargin<3
    options = [];
end

if strcmpi(Type,'MST')
    DistMatrix = sparse(DistMatrix);
    DistMatrix = tril(DistMatrix,-1);
    [ST, PRED] = graphminspantree(DistMatrix,'Method','Kruskal');
elseif strcmpi(Type,'LAST')
    DistMatrix = sparse(DistMatrix);
    TrilDistMatrix = tril(DistMatrix, -1);
    alpha = getoptions(options,'alpha',8);
    if ~isnumeric(alpha)
        [ST,PRED] = graphminspantree(TrilDistMatrix,'Method','Kruskal');
        RootNode = find(PRED==0);
        [DistRootJ,~,~] = graphshortestpath(ST,RootNode,'directed',false);
        PointwiseDistortion = DistRootJ./DistMatrix(RootNode,:);
        alpha = mean(PointwiseDistortion);
    end
    
    OPT_WEIGHT = Inf;
    RowNodes = 1:size(DistMatrix,1);
    Weights = zeros(size(RowNodes));
    for RootNode=1:size(DistMatrix,1)
        ST = graphminspantree(TrilDistMatrix,RootNode,'Method','Kruskal');
        [DistRootJ,~,~] = graphshortestpath(ST,RootNode,'directed',false);
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
    PRED = OPT_PRED;
    RowNodes(OPT_RootNode) = [];
    OPT_PRED(OPT_RootNode) = [];
    OPT_Weights(OPT_RootNode) = [];
    SPT = sparse(RowNodes,OPT_PRED,OPT_Weights,size(DistMatrix,1),size(DistMatrix,2));
    ST = max(SPT,SPT');
    ST = tril(ST, -1);
else
    error('Unidefined Graph Type!');
end

end

