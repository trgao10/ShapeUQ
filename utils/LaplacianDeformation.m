function [TextureCoords1,ImprMap] = LaplacianDeformation(GM,GN,InputMap,FeatureType,TextureCoords1,TextureCoords2,options)
%LAPLACIANDEFORMATION: Laplacian mesh deformation
%   Current design prefers geodesic mutually nearest neighboring features
%   than Euclidean mutually nearest neighboring features: if a feature is
%   assigned to different target features under these two ad-hoc
%   assignments, then the "Euclidean target" gets "wiped out" in the
%   "unique" section.

if nargin<7
    options = [];
end
TextureCoords2_kdtree = getoptions(options,'TextureCoords2_kdtree',kdtree_build(TextureCoords2'));
GaussMinMatch = getoptions(options,'GaussMinMatch','on');
debug = getoptions(options,'debug','off');

%%% extract matching features
[~,LAP_FEATURESN,preLAP_FEATURESM] = FindMutuallyNearestNeighbors(GM,GN,InputMap,FeatureType);
[~,LAP_EUC_FEATURESN,preLAP_EUC_FEATURESM] = FindEuclideanMutuallyNearestNeighbors(GM,GN,InputMap,FeatureType);
LAP_FEATURESM_INDS = [preLAP_FEATURESM;preLAP_EUC_FEATURESM];
LAP_FEATURESN_INDS = [LAP_FEATURESN;LAP_EUC_FEATURESN];
[LAP_FEATURESM_INDS,NoRepeatInds] = unique(LAP_FEATURESM_INDS);
LAP_FEATURESN_INDS = LAP_FEATURESN_INDS(NoRepeatInds);
if strcmpi(GaussMinMatch,'on') && ~strcmpi(FeatureType,'GaussMin')
    [~,LAP_GAUSSMIN_FEATURESN,preLAP_GAUSSMIN_FEATURESM] = FindMutuallyNearestNeighbors(GM,GN,InputMap,'GaussMin');
    LAP_FEATURESM_INDS = [LAP_FEATURESM_INDS;preLAP_GAUSSMIN_FEATURESM];
    LAP_FEATURESN_INDS = [LAP_FEATURESN_INDS;LAP_GAUSSMIN_FEATURESN];
end

%%% Laplacian deformation
% BV = GM.FindBoundaries;
% L = GM.Aux.LB;
%%%%% use normalized combinatorial laplacian? 
% D = sum(GM.A);
% L = speye(GM.nV)-sparse(diag(1./D))*GM.A;
%%%%% use planar contangent laplacian?
if strcmpi(debug,'on')
    figure;subplot(1,2,1);
    TR = triangulation(GM.F',TextureCoords1');
    triplot(TR);axis equal;hold on;
    scatter(TextureCoords1(1,LAP_FEATURESM_INDS),TextureCoords1(2,LAP_FEATURESM_INDS),30,'r','filled');
    scatter(TextureCoords2(1,LAP_FEATURESN_INDS),TextureCoords2(2,LAP_FEATURESN_INDS),30,'k','filled');
end

KM = Mesh('VF',TextureCoords1,GM.F);
BV = KM.FindBoundaries;
L = KM.ComputeCotanLaplacian;
L([LAP_FEATURESM_INDS;BV],:) = 0;
L([LAP_FEATURESM_INDS;BV]+([LAP_FEATURESM_INDS;BV]-1)*size(L,1)) = 1;
rhs = zeros(size(TextureCoords1,2),2);
rhs(BV,:) = TextureCoords1(:,BV)';
rhs(LAP_FEATURESM_INDS,:) = TextureCoords2(:,LAP_FEATURESN_INDS)';
TextureCoords1 = (L\rhs)';

if strcmpi(debug,'on')
    subplot(1,2,2);
    TR = triangulation(GM.F',TextureCoords1');
    triplot(TR);axis equal;hold on;
    scatter(TextureCoords1(1,LAP_FEATURESM_INDS),TextureCoords1(2,LAP_FEATURESM_INDS),30,'r','filled');
    scatter(TextureCoords2(1,LAP_FEATURESN_INDS),TextureCoords2(2,LAP_FEATURESN_INDS),30,'k','filled');
    keyboard
end

%%% construct vertex permutation map
ImprMap = kdtree_nearest_neighbor(TextureCoords2_kdtree, TextureCoords1');

end

