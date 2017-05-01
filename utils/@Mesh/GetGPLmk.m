function [GPLmkIdx,ptuq] = GetGPLmk(G,numLmk)
%GETGPLMK Summary of this function goes here
%   Detailed explanation goes here

% if nargin < 3
%     lambda = 0.5;
% end

G.Centralize('ScaleArea');
[~,TriArea] = G.ComputeSurfaceArea();
G.Aux.VertArea = G.F2V'*TriArea;

% [Cgauss,Cmean] = G.ComputeCurvature();
% Lambda = G.Aux.VertArea.*(lambda*abs(Cgauss)/sum(abs(Cgauss))+(1-lambda)*abs(Cmean)/sum(abs(Cmean)));

[~,curvature] = findPointNormals(G.V',10);
Lambda = G.Aux.VertArea.*curvature/sum(curvature);

if size(G.E,1) > 2
    [I,J] = find(tril(G.E));
    G.E = ([I,J])';
end
EdgeIdxI = G.E(1,:);
EdgeIdxJ = G.E(2,:);
bandwidth = mean(sqrt(sum((G.V(:,EdgeIdxI)-G.V(:,EdgeIdxJ)).^2)))/3;

BNN = min(500,G.nV);
atria = nn_prepare(G.V');
[idx, dist] = nn_search(G.V',atria,(1:G.nV)',BNN+1,-1,0.0);
fullPhi = sparse(repmat(1:G.nV,1,BNN+1),idx,exp(-dist.^2/bandwidth),G.nV,G.nV);

% PDistMat = squareform(pdist(G.V'));
% fullPhi = exp(-PDistMat.^2/bandwidth);

disp('Constructing full kernel......');
tic;
fullMatProd = fullPhi * sparse(1:G.nV,1:G.nV,Lambda,G.nV,G.nV) * fullPhi;
disp(['full kernel constructed in ' num2str(toc) ' sec.']);

KernelTrace = diag(fullMatProd);
GPLmkIdx = zeros(1,numLmk);

cback = 0;
for j=1:numLmk
    for cc=1:cback
        fprintf('\b');
    end
    cback = fprintf('Landmark: %4d\n',j);
    
    if j == 1
        ptuq = KernelTrace;
    else
        ptuq = KernelTrace - sum(fullMatProd(:,GPLmkIdx(1:(j-1)))'...
            .*(fullMatProd(GPLmkIdx(1:(j-1)),GPLmkIdx(1:(j-1)))\fullMatProd(GPLmkIdx(1:(j-1)),:)),1)';
    end
    [~,maxUQIdx] = max(ptuq);
    GPLmkIdx(j) = maxUQIdx;
end

ptuq = KernelTrace - sum(fullMatProd(GPLmkIdx,:)...
    .*(fullMatProd(GPLmkIdx,GPLmkIdx)\fullMatProd(GPLmkIdx,:)),1)';

end

