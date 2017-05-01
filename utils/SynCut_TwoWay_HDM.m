function [rslt] = SynCut_TwoWay_HDM(G, edgePotCell, params)
%SYNCUT_TWOWAY_HDM
%   INPUTS:
%     G ---------------------------- graph
%     edgePotCell ------------------ edge potential in cell array format
%     params ----------------------- struct holding all parameters
%
%   OUTPUTS:
%     rslt ------------------------- struct holding all outputs
%
%
% Tingran Gao (trgao10@math.duke.edu)
% last modified: Oct 16, 2016
%

debugFlag = getoptions(params, 'debugFlag', false);
d = getoptions(params, 'd', Inf);
numClusters = getoptions(params, 'numClusters', 2);
tol = getoptions(params, 'tol', 1e-8);
maxIter = getoptions(params, 'maxIter', 10);
numKmeans = getoptions(params, 'numKmeans', 200);
bandwidth = getoptions(params, 'bandwidth', 1);
adjType = getoptions(params, 'adjType', 'dis');

if numClusters > 2
    error('This routine only demos the usage for the case numCluster=2');
end

[GCL, GCL_Dvec, GCL_W] = assembleGCL(G.adjMat, edgePotCell, d);

if debugFlag
    hsv = getoptions(params, 'hsv', rgb2hsv(winter));
    close(gcf);
    try
        vertPotCell = params.vertPotCell;
    catch
        error('with debugFlag=true, must provide ground truth vertPotCell');
    end
    try
        colorList = getoptions(params, 'colorList', {'r','b','k','m'});
    catch
        error('with debugFlag=true, must provide ground truth colorList');
    end
    
    figure('Position',[30,550,560,420]);
    plot(graph(G.adjMat), 'XData', G.V(:,1), 'YData', G.V(:,2), 'LineWidth', 0.1);
    axis equal
    axis([0,numClusters,0,1]);
    hold on
    for j=1:length(G.ccRowIdx)
        staPtCoords = G.V(G.ccRowIdx(j),:);
        endPtCoords = G.V(G.ccColIdx(j),:);
        line([staPtCoords(:,1);endPtCoords(:,1)],...
            [staPtCoords(:,2);endPtCoords(:,2)],'Color','r','LineStyle','-');
    end
    title(sprintf('Spectral Gap = %.2f, CCE/TTE = %d/%d',...
        G.specGap, length(G.ccRowIdx), nnz(G.adjMat)/2),'Interpreter','latex');

    vertPotMat = cat(1, vertPotCell{:});
    GCL_unnormalized = diag(GCL_Dvec) - GCL_W;
    [GroundTruthPerEdgeFrustVec, GroundTruthPerEdgeFrustMat] =...
        getPerEdgeFrustFromEdgePot(G.adjMat, edgePotCell, vertPotCell);
    rescaled_vertPotMat = diag(sqrt(GCL_Dvec))*vertPotMat;
    fprintf('[GroundTruth] Rayleigh quotient (normalized GCL) = %f\n',...
        trace(rescaled_vertPotMat'*GCL*rescaled_vertPotMat));
    fprintf('[GroundTruth] Rayleigh quotient (unnormalized GCL) = %f\n',...
        trace(vertPotMat'*GCL_unnormalized*vertPotMat));
    fprintf('[GroundTruth] total frustration = %f\n',...
        sum(GroundTruthPerEdgeFrustVec));
end

xi = Inf;
iterCounter = 0;
wAdjMat = G.adjMat;
while true
    iterCounter = iterCounter+1;
    if debugFlag
        fprintf('++++++++ Iteration %d +++++++++++++++++++++++++++++\n',...
                iterCounter);
    end
    
    xi_old = xi;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Step 1. Synchronization (global)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [wGCL, wGCL_Dvec, wGCL_W] = assembleGCL(wAdjMat, edgePotCell, d);
    [RelaxSolCell, RelaxSolMat] = syncSpecRelax(wGCL, d, wGCL_Dvec);
    [RelaxSolPerEdgeFrustVec, RelaxSolPerEdgeFrustMat] =...
        getPerEdgeFrustFromEdgePot(G.adjMat, edgePotCell, RelaxSolCell);
    if iterCounter == 1
        origRelaxSolCell = RelaxSolCell;
    end
    
    if debugFlag
        rescaled_RelaxSolMat = diag(sqrt(GCL_Dvec))*RelaxSolMat;
        fprintf('[RelaxSol] Rayleigh quotient (normalized GCL) = %f\n',...
            trace(rescaled_RelaxSolMat'*GCL*rescaled_RelaxSolMat));
        fprintf('[RelaxSol] Rayleigh quotient (unnormalized GCL) = %f\n',...
            trace(RelaxSolMat'*GCL_unnormalized*RelaxSolMat));
        fprintf('[RelaxSol] total frustration = %f\n',...
            sum(RelaxSolPerEdgeFrustVec));        
        figure;
        subplot(1,2,1);
        imagesc(RelaxSolPerEdgeFrustMat);
        axis square
        title(sprintf('RelaxSolPerEdgeFrustMat: %.2f', sum(RelaxSolPerEdgeFrustVec)));
        subplot(1,2,2);
        imagesc(G.adjMat);
        axis square
        title('adjacency matrix');
%         keyboard
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Step 2. Spectral Clustering, Pass 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while true
        clusterLabel =...
            specClusterWrapper(RelaxSolPerEdgeFrustMat,numClusters,...
            numKmeans,bandwidth,adjType);
        %%%%% we expect spectral clustering to partition the graph into two
        %%%%% connected components
        connTest = zeros(size(clusterLabel));
        numConnComp = zeros(size(clusterLabel));
        for j=1:length(connTest)
            connTest(j) = ~isempty(find(sum(G.adjMat(clusterLabel{j},clusterLabel{j})) == 0, 1));
            numConnComp(j) = length(unique(conncomp(graph(G.adjMat(clusterLabel{j},clusterLabel{j})))));
        end
        if (sum(connTest) == 0) && (all(numConnComp == 1))
            break
        end
        error('error: spectral clustering does not produce connected components');
    end
                       
    if debugFlag
        %%%% for more consistent visualization effect, always set left
        %%%% cluster to be red an right cluster to be blue
        %%%% Note this will only work for K=2 (binary clustering)!
        n = size(G.adjMat,1) / numClusters;
        if sum(clusterLabel{1} <= n) < sum(clusterLabel{2} <= n)
            tmp = clusterLabel{1};
            clusterLabel{1} = clusterLabel{2};
            clusterLabel{2} = tmp;
            clear tmp
        end
        
        %%%% plot spectral clustering results
        if ~exist('specClusteringFigure', 'var')
            specClusteringFigure = figure('Position',[900,150,1000,400]);
        else
            figure(specClusteringFigure);
        end
        subplot(1,2,1);
        for j=1:length(clusterLabel)
            scatter(G.V(clusterLabel{j},1),G.V(clusterLabel{j},2),...
                    20,colorList{j},'filled');
            if j==1
                hold on
            end
        end
        axis equal
        axis([0,numClusters,0,1]);
        line([1,1],[0,1],'Color','g');
        title(sprintf('Spectral Clustering in Iteration %d, Pass 1',iterCounter));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Step 3. Synchronization (local)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rescaled_relaxSol_cluster = cell(1,numClusters);
%     try
        for j=1:numClusters
            adjMat_cluster = G.adjMat(clusterLabel{j},clusterLabel{j});
            edgePotCell_cluster = edgePotCell(clusterLabel{j},clusterLabel{j});
            [GCL_cluster, GCL_Dvec_cluster, GCL_W_cluster] = assembleGCL(adjMat_cluster, edgePotCell_cluster, d);
            [~, rescaled_relaxSol_cluster{j}] = syncSpecRelax(GCL_cluster, d, GCL_Dvec_cluster);
            
            if debugFlag
                GCL_unnormalized_cluster = diag(GCL_Dvec_cluster) - GCL_W_cluster;
                fprintf('[cluster %d] rescaled relaxed solution Rayleigh quotient (unnormalized Laplacian) = %f\n',...
                    j, trace(rescaled_relaxSol_cluster{j}'*GCL_unnormalized_cluster*rescaled_relaxSol_cluster{j}));
            end
        end
%     catch ME
%         disp(ME.message);
%         keyboard
%         return
%     end
    
    combinedSolCell = cell(1,size(G.adjMat,1));
    combinedSolMat = zeros(size(G.adjMat,1)*d,d);
    for j=1:numClusters
        for k=1:length(clusterLabel{j})
            tmpIdx = ((clusterLabel{j}(k)-1)*d+1):(clusterLabel{j}(k)*d);
            tmpIdxInCluster = ((k-1)*d+1):(k*d);
            tmpBlock = rescaled_relaxSol_cluster{j}(tmpIdxInCluster,:);
            [U,S,V] = svd(tmpBlock);
            combinedSolCell{clusterLabel{j}(k)} = U*V';
            combinedSolMat(tmpIdx,:) = combinedSolCell{clusterLabel{j}(k)};
        end
    end
    
    [combinedSolPerEdgeFrustVec, combinedSolPerEdgeFrustMat] =...
            getPerEdgeFrustFromEdgePot(G.adjMat, edgePotCell, combinedSolCell);
%     if debugFlag        
%         rescaled_combinedSolMat = diag(sqrt(GCL_Dvec))*combinedSolMat;
%         fprintf('[combinedSol] Rayleigh quotient (normalized Laplacian) = %f\n',...
%             trace(rescaled_combinedSolMat'*GCL*rescaled_combinedSolMat));
%         fprintf('[combinedSol] Rayleigh quotient (unnormalized Laplacian) = %f\n',...
%             trace(combinedSolMat'*GCL_unnormalized*combinedSolMat));
%         fprintf('[combinedSol] total frustration = %f\n',...
%             sum(combinedSolPerEdgeFrustVec));
%     end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Step 4. Collage
    % [ATTENTION] Here we only deal with the simplest case of 2 clusters!!
    % [TODO] Try implementing minimum instead of average?
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CollageSolCell = combinedSolCell;
    CollageSolMat = combinedSolMat;
    
    %%%%%%%%% extract cross-cluster edges
    crossClusterAdjMat = G.adjMat;
    for j=1:numClusters
        crossClusterAdjMat(clusterLabel{j},clusterLabel{j}) = 0;
    end
    [ccRowIdx, ccColIdx] = find(triu(crossClusterAdjMat));
    
    %%%%%%%%% orient all cross-cluster edges to go from lower to higher
    for j=1:length(ccRowIdx)
        if isempty(find(clusterLabel{1} == ccRowIdx(j), 1))
            tmpIdx = ccRowIdx(j);
            ccRowIdx(j) = ccColIdx(j);
            ccColIdx(j) = tmpIdx;
        end
    end
    
    %%%%%%%%% check all cross-cluster edges re-oriented
    counter = 0;
    for j=1:length(ccRowIdx)
        if isempty(find(clusterLabel{1} == ccRowIdx(j), 1))
            counter = counter + 1;
        end
    end
    if counter > 0
        fprintf('counter = %d\n', counter);
        error('edge oriented reversely');
    end
    
    %%%%%%%%% form ensemble covariance matrix and SVD
    ccCorr = zeros(d);
    for j=1:length(ccRowIdx)
        blockRowIdx = ((ccRowIdx(j)-1)*d+1):(ccRowIdx(j)*d);
        blockColIdx = ((ccColIdx(j)-1)*d+1):(ccColIdx(j)*d);
%         ccCorr = ccCorr+RelaxSolPerEdgeFrustMat(ccRowIdx(j),ccColIdx(j))*CollageSolCell{ccColIdx(j)}'*GCL_W(blockRowIdx, blockColIdx)'*CollageSolCell{ccRowIdx(j)};
        ccCorr = ccCorr+combinedSolPerEdgeFrustMat(ccRowIdx(j),ccColIdx(j))*CollageSolCell{ccColIdx(j)}'*GCL_W(blockRowIdx, blockColIdx)'*CollageSolCell{ccRowIdx(j)};
    end
    [U,~,V] = svd(ccCorr);
    collageR = V*U';
    
    %%%%%%%%% perform the collage
    for k=1:length(clusterLabel{1})
        CollageSolCell{clusterLabel{1}(k)} = CollageSolCell{clusterLabel{1}(k)}*collageR;
        tmpIdx = ((clusterLabel{1}(k)-1)*d+1):(clusterLabel{1}(k)*d);
        CollageSolMat(tmpIdx,:) = CollageSolMat(tmpIdx,:)*collageR;
    end
    
    [CollageSolPerEdgeFrustVec, CollageSolPerEdgeFrustMat] =...
                       getPerEdgeFrustFromEdgePot(G.adjMat, edgePotCell, CollageSolCell);
    if debugFlag
        inClusterTotalFrust = zeros(1,2);
        bwClusterTotalFrust = zeros(1,2);
        bwClusterTotalFrustVecs = zeros(2,length(ccRowIdx));
        for k=1:length(ccRowIdx)
            bwClusterTotalFrust(1) = bwClusterTotalFrust(1)+combinedSolPerEdgeFrustMat(ccRowIdx(k),ccColIdx(k));
            bwClusterTotalFrustVecs(1,k) = combinedSolPerEdgeFrustMat(ccRowIdx(k),ccColIdx(k));
            bwClusterTotalFrust(2) = bwClusterTotalFrust(2)+CollageSolPerEdgeFrustMat(ccRowIdx(k),ccColIdx(k));
            bwClusterTotalFrustVecs(2,k) = CollageSolPerEdgeFrustMat(ccRowIdx(k),ccColIdx(k));
        end
        inClusterTotalFrustCells = cell(2,length(clusterLabel));
        for j=1:length(clusterLabel)
            inClusterTotalFrust(1) = inClusterTotalFrust(1)+sum(sum(combinedSolPerEdgeFrustMat(clusterLabel{j},clusterLabel{j})))/2;
            [~,~,inClusterTotalFrustCells{1,j}] = find(triu(combinedSolPerEdgeFrustMat(clusterLabel{j},clusterLabel{j})));
            inClusterTotalFrust(2) = inClusterTotalFrust(2)+sum(sum(CollageSolPerEdgeFrustMat(clusterLabel{j},clusterLabel{j})))/2;
            [~,~,inClusterTotalFrustCells{2,j}] = find(triu(CollageSolPerEdgeFrustMat(clusterLabel{j},clusterLabel{j})));
        end
        
        rescaled_combinedSolMat = diag(sqrt(GCL_Dvec))*combinedSolMat;
        fprintf('[combinedSol] Rayleigh quotient (normalized Laplacian) = %f\n',...
            trace(rescaled_combinedSolMat'*GCL*rescaled_combinedSolMat));
        fprintf('[combinedSol] total frustration = %f\n',...
            sum(combinedSolPerEdgeFrustVec));
        fprintf('[combinedSol] (in+bw clusters) frustration = %f\n',...
            sum(inClusterTotalFrust(1)+bwClusterTotalFrust(1)));
        fprintf('[combinedSol] in-cluster frustration = %f\n',...
            sum(inClusterTotalFrust(1)));
        fprintf('[combinedSol] bw-cluster frustration = %f\n',...
            sum(bwClusterTotalFrust(1)));
        
        rescaled_CollageSolMat = diag(sqrt(GCL_Dvec))*CollageSolMat;
        fprintf('[CollageSol] Rayleigh quotient (normalized Laplacian) = %f\n',...
            trace(rescaled_CollageSolMat'*GCL*rescaled_CollageSolMat));
        fprintf('[CollageSol] total frustration = %f\n',...
            sum(CollageSolPerEdgeFrustVec));
        fprintf('[CollageSol] (in+bw clusters) frustration = %f\n',...
            sum(inClusterTotalFrust(2)+bwClusterTotalFrust(2)));
        fprintf('[CollageSol] in-cluster frustration = %f\n',...
            sum(inClusterTotalFrust(2)));
        fprintf('[CollageSol] bw-cluster frustration = %f\n',...
            sum(bwClusterTotalFrust(2)));
        
        if ~exist('inCluserFigure', 'var')
            inCluserFigure = figure('Name', 'in-cluster edge-wise frustrations should coincide',...
                'Position', [1150,500,560,420]);
        else
            figure(inCluserFigure);
        end
        subplot(1,2,1);
        hist(cat(1,inClusterTotalFrustCells{1,:}));
        title(sprintf('combinedSol, %.2f',sum(cat(1,inClusterTotalFrustCells{1,:}))));
        subplot(1,2,2);
        hist(cat(1,inClusterTotalFrustCells{2,:}));
        title(sprintf('CollageSol, %.2f',sum(cat(1,inClusterTotalFrustCells{2,:}))));
        
        if ~exist('bwClusterFigure', 'var')
            bwClusterFigure = figure('Name', 'bw-cluster edge-wise frustrations',...
                'Position', [1300,500,560,420]);
        else
            figure(bwClusterFigure);
        end
        subplot(1,2,1);
        hist(bwClusterTotalFrustVecs(1,:));
        title(sprintf('combinedSol, %.2f',sum(bwClusterTotalFrustVecs(1,:))));
        subplot(1,2,2);
        hist(bwClusterTotalFrustVecs(2,:));
        title(sprintf('CollageSol, %.2f',sum(bwClusterTotalFrustVecs(2,:))));
        
        if ~exist('perEdgeFrustFigure', 'var')
            perEdgeFrustFigure = figure('Position',[50,100,800,600]);
        else
            figure(perEdgeFrustFigure);
        end
        subplot(2,2,1);
        plotPerEdgeFrustration_enhance_cc(G,GroundTruthPerEdgeFrustMat,hsv);
        title(['GroundTruth total frustration = '...
            num2str(sum(GroundTruthPerEdgeFrustVec))]);
        subplot(2,2,2);
        plotPerEdgeFrustration_enhance_cc(G,RelaxSolPerEdgeFrustMat,hsv);
        title(['RelaxSol total frustration = '...
            num2str(sum(RelaxSolPerEdgeFrustVec))]);
        subplot(2,2,3);
        plotPerEdgeFrustration_enhance_cc(G,combinedSolPerEdgeFrustMat,hsv);
        title(['combinedSol total frustration = '...
            num2str(sum(combinedSolPerEdgeFrustVec))]);
        subplot(2,2,4);
        plotPerEdgeFrustration_enhance_cc(G,CollageSolPerEdgeFrustMat,hsv);
        title(['CollageSol total frustration = '...
            num2str(sum(CollageSolPerEdgeFrustVec))]);
        pause();
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Step 5. Spectral Clustering, Pass 2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clusterLabel =...
        specClusterWrapper(CollageSolPerEdgeFrustMat,numClusters,...
                           numKmeans,bandwidth,adjType);
    xi = evalXi(G, clusterLabel, d, CollageSolPerEdgeFrustMat);
    
    if debugFlag
        %%%% for more consistent visualization effect, always set left
        %%%% cluster to be red an right cluster to be blue
        %%%% Note this will only work for K=2 (binary clustering)!
        n = size(G.adjMat,1) / numClusters;
        if sum(clusterLabel{1} <= n) < sum(clusterLabel{2} <= n)
            tmp = clusterLabel{1};
            clusterLabel{1} = clusterLabel{2};
            clusterLabel{2} = tmp;
            clear tmp
        end
        
        %%%% plot spectral clustering results
        figure(specClusteringFigure);
        subplot(1,2,2);
        for j=1:length(clusterLabel)
            scatter(G.V(clusterLabel{j},1),G.V(clusterLabel{j},2),...
                    20,colorList{j},'filled');
            if j==1
                hold on
            end
        end
        axis equal
        axis([0,numClusters,0,1]);
        line([1,1],[0,1],'Color','g');
        title(sprintf('Spectral Clustering in Iteration %d, Pass 2',iterCounter));
        
        fprintf('iteration %d, xi = %4d\n', iterCounter, xi);
        fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++\n');
        pause();
    end
    
    [wRowIdx,wColIdx,wVals] = find(CollageSolPerEdgeFrustMat);
    wAdjMat = sparse(wRowIdx,wColIdx,exp(-wVals/mean(wVals)),size(wAdjMat,1),size(wAdjMat,2));
    
    if (abs(xi) < tol) || (iterCounter >= maxIter) || ((xi_old < Inf) && (abs(xi-xi_old) < tol*xi_old))
        break
    end
end

%%%%%% very ad-hoc error counts and error rates computation
nPtsPerCluster = size(G.adjMat, 1)/numClusters;
errCounts = min(sum(clusterLabel{1} > nPtsPerCluster)+sum(clusterLabel{2} <= nPtsPerCluster),...
                sum(clusterLabel{1} <= nPtsPerCluster)+sum(clusterLabel{2} > nPtsPerCluster));
errRate = errCounts / size(G.adjMat, 1);

rslt = struct('G', G, 'params', params, 'iterCounter', iterCounter,...
              'errCounts', errCounts, 'errRate', errRate);
[rslt.edgePotCell] = deal(edgePotCell);
[rslt.clusterLabel] = deal(clusterLabel);
[rslt.CollageSolCell] = deal(CollageSolCell);
[rslt.RelaxSolCell] = deal(origRelaxSolCell);

if exist('combinedSolCell', 'var')
    [rslt.combinedSolCell] = deal(combinedSolCell);
end

end
