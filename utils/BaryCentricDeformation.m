function [TextureCoords1,ImprMap] = BaryCentricDeformation(GM,GN,InputMap,FeatureType,TextureCoords1,TextureCoords2,options)
%BARYCENTRICDEFORMATION Summary of this function goes here
%   Detailed explanation goes here

if nargin<7
    options = [];
end
TextureCoords2_kdtree = getoptions(options,'TextureCoords2_kdtree',kdtree_build(TextureCoords2'));
GaussMinMatch = getoptions(options,'GaussMinMatch','on');

%%% extract matching features
[~,TPS_FEATURESN,preTPS_FEATURESM] = FindMutuallyNearestNeighbors(GM,GN,InputMap,FeatureType);
[~,TPS_EUC_FEATURESN,preTPS_EUC_FEATURESM] = FindEuclideanMutuallyNearestNeighbors(GM,GN,InputMap,FeatureType);
TPS_FEATURESM_INDS = [preTPS_FEATURESM;preTPS_EUC_FEATURESM];
TPS_FEATURESN_INDS = [TPS_FEATURESN;TPS_EUC_FEATURESN];
[TPS_FEATURESM_INDS,NoRepeatInds] = unique(TPS_FEATURESM_INDS);
TPS_FEATURESN_INDS = TPS_FEATURESN_INDS(NoRepeatInds);
if strcmpi(GaussMinMatch,'on') && ~strcmpi(FeatureType,'GaussMin')
    [~,TPS_GAUSSMIN_FEATURESN,preTPS_GAUSSMIN_FEATURESM] = FindMutuallyNearestNeighbors(GM,GN,InputMap,'GaussMin');
    TPS_FEATURESM_INDS = [TPS_FEATURESM_INDS;preTPS_GAUSSMIN_FEATURESM];
    TPS_FEATURESN_INDS = [TPS_FEATURESN_INDS;TPS_GAUSSMIN_FEATURESN];
end

BV = GM.FindBoundaries;
anchors = TextureCoords1(:,[BV;TPS_FEATURESM_INDS]);
target_anchors = [TextureCoords1(:,BV),TextureCoords2(:,TPS_FEATURESN_INDS)];
Tri = delaunayTriangulation(anchors');
NewTextureCoords1 = zeros(size(TextureCoords1));

nV = size(NewTextureCoords1,2);
nF = size(Tri.ConnectivityList,1);
for j = 1:nV
    progressbar(j,nV,20);
    BC = Tri.cartesianToBarycentric((1:nF)',repmat(TextureCoords1(:,j)',nF,1));
    tind = find(all(BC>-1e-10,2));
    if(numel(tind)>1)
        tind = tind(1);
    elseif numel(tind)==0
        warning(['Point ', num2str(j), ' was not found in any triangle. It stays where it was.']);
        NewTextureCoords1(:,j) = TextureCoords1(:,j);
        continue;
    end
    BC = BC(tind,:);
    NewTextureCoords1(:,j) = target_anchors(:,Tri.ConnectivityList(tind,:))*BC';
end
TextureCoords1 = NewTextureCoords1;

%%% construct vertex permutation map
ImprMap = kdtree_nearest_neighbor(TextureCoords2_kdtree, TextureCoords1');

end

