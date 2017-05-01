function [rslt] = FeatureFix(GM,GN,TAXAind1,TAXAind2,options)
%FEATUREFIX:  fix features after composing maps along a propagation path
%   rslt.Gname1:            name of the first mesh
%   rslt.Gname2:            name of the second mesh
%   rslt.ImprDist:            continuous Procrustes distance
%   rslt.ImprMap:             optimal map generating cP distance
%   rslt.invImprMap:          inverse of rslt.ImprMap
%   rslt.TextureCoords1:    texture coordinates for the first mesh
%                           (deformed)
%   rslt.TextureCoords2:    textrue coordinates for the second mesh
%                           (not deformed)
%   rslt.ref:               =0 if Improved map is orientation-preserving
%                           =1 if Improved map is orientation-reversing
%
%   Tingran Gao, trgao10@math.duke.edu
%   last modified: 23 Aug 2014
%

if nargin<5
    options = [];
end
ChunkSize = getoptions(options,'ChunkSize',55);
GroupSize = getoptions(options,'GroupSize',116);
if ~isfield(GM.Aux,'name') && ~isfield(GN.Aux,'name')
    error('Either Mesh missing .Aux.name!');
end
if ~isfield(options,'TextureCoords1Path') && ~isfield(options,'TextureCoords2Path')
    error('TextureCoords1Path or TextureCoords1Path missing!');
end

ChunkIdx = @(TAXAind1,TAXAind2) ceil(((TAXAind1-1)*GroupSize+TAXAind2)/ChunkSize);

rslt.Gname1 = GM.Aux.name;
rslt.Gname2 = GN.Aux.name;

load([options.TextureCoords1Path 'TextureCoords1_mat_' num2str(ChunkIdx(TAXAind1,TAXAind2)) '.mat']);
load([options.TextureCoords2Path 'TextureCoords2_mat_' num2str(ChunkIdx(TAXAind1,TAXAind2)) '.mat']);
rslt.TextureCoords1 = TextureCoords1Matrix{TAXAind1,TAXAind2};
rslt.TextureCoords2 = TextureCoords2Matrix{TAXAind1,TAXAind2};

MapMN = knnsearch(rslt.TextureCoords2',rslt.TextureCoords1');
disp('Performing Feature Fixing...');
[rslt.TextureCoords1,rslt.ImprMap] = TPSDeformation(GM,GN,MapMN,'ConfMax',rslt.TextureCoords1,rslt.TextureCoords2);
disp('Done.');

rslt.ImprMap = knnsearch(rslt.TextureCoords2',rslt.TextureCoords1');
[rslt.ImprDist,R] = MapToDist(GM.V,GN.V,rslt.ImprMap,GM.Aux.VertArea);
if det(R)>0
    rslt.ref = 0;
else
    rslt.ref = 1;
end
rslt.invImprMap = knnsearch(rslt.TextureCoords1',rslt.TextureCoords2');

end

