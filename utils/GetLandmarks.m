function [LandmarkInds,Landmarks] = GetLandmarks(Name,LandmarksPath,MeshPath,options)
%GETLANDMARKS Summary of this function goes here
%   Detailed explanation goes here

if nargin<4
    options = [];
end
NumLandmark = getoptions(options,'NumLandmark',16);

G = Mesh('off',MeshPath);

LandmarkFile = load(LandmarksPath);
rawLandmarks = LandmarkFile.PP(strcmpi(LandmarkFile.names,Name),1:NumLandmark,:);
Landmarks = zeros(size(rawLandmarks,2),3);
for k=1:size(rawLandmarks,2)
    Landmarks(k,:) = [rawLandmarks(1,k,1), rawLandmarks(1,k,2), rawLandmarks(1,k,3)];
end
if ~isempty(strfind(LandmarksPath,'Clement'))
    tree = KDTreeSearcher(G.V');
    LandmarkInds = tree.knnsearch(Landmarks);
    Landmarks = G.V(:,LandmarkInds)';
%     Landmarks = (Landmarks-repmat(G.Aux.Center',NumLandmark,1))*sqrt(1/G.Aux.Area);
else
    [G.Aux.Area,G.Aux.Center] = G.Centralize('ScaleArea');
    tree = KDTreeSearcher(G.V');
    Landmarks = Landmarks-repmat(G.Aux.Center',NumLandmark,1)*sqrt(1/G.Aux.Area);
    LandmarkInds = tree.knnsearch(Landmarks);
    Landmarks = G.V(:,LandmarkInds)';
end
% Landmarks = (Landmarks*sqrt(G.Aux.Area)-repmat(G.Aux.Center',NumLandmark,1))*sqrt(1/G.Aux.Area);

end

