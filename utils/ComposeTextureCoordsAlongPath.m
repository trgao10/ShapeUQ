function [TextureCoords1,TextureCoords2] = ComposeTextureCoordsAlongPath(InputPath,TextureCoords1Path,TextureCoords2Path,GroupSize,ChunkSize,ProgressBar)
%COMPOSETEXTURECOORDSALONGPATH Summary of this function goes here
%   Detailed explanation goes here

ChunkIdx = @(TAXAind1,TAXAind2) ceil(((TAXAind1-1)*GroupSize+TAXAind2)/ChunkSize);

if length(InputPath)==1
    load([TextureCoords1Path 'TextureCoords1_mat_' num2str(ChunkIdx(InputPath(1),InputPath(1))) '.mat']);
    load([TextureCoords2Path 'TextureCoords2_mat_' num2str(ChunkIdx(InputPath(1),InputPath(1))) '.mat']);
    TextureCoords1 = TextureCoords1Matrix{InputPath(1),InputPath(1)};
    TextureCoords2 = TextureCoords2Matrix{InputPath(1),InputPath(1)};
    return;
end

load([TextureCoords1Path 'TextureCoords1_mat_' num2str(ChunkIdx(InputPath(end-1),InputPath(end))) '.mat']);
load([TextureCoords2Path 'TextureCoords2_mat_' num2str(ChunkIdx(InputPath(end-1),InputPath(end))) '.mat']);
TextureCoords1 = TextureCoords1Matrix{InputPath(end-1),InputPath(end)};
TextureCoords2 = TextureCoords2Matrix{InputPath(end-1),InputPath(end)};

if (length(InputPath)>2)
    for j=3:length(InputPath)
        load([TextureCoords1Path 'TextureCoords1_mat_' num2str(ChunkIdx(InputPath(end-j+1),InputPath(end-j+2))) '.mat']);
        load([TextureCoords2Path 'TextureCoords2_mat_' num2str(ChunkIdx(InputPath(end-j+1),InputPath(end-j+2))) '.mat']);
        NextNodeTextureCoords1 = TextureCoords1Matrix{InputPath(end-j+1),InputPath(end-j+2)};
        NextNodeTextureCoords2 = TextureCoords2Matrix{InputPath(end-j+1),InputPath(end-j+2)};
        TextureCoords1 = PropagateTextureCoords(TextureCoords1,NextNodeTextureCoords1,NextNodeTextureCoords2,ProgressBar);
    end
end

end

function NewTextureCoords1 = PropagateTextureCoords(TextureCoords1,NextTextureCoords1,NextTextureCoords2,ProgressBar)
% NextTextureCoords2 is the un-deformed version of TextureCoords1.
% This function returns NextTextureCoords1 under the same deformation that
% deformed NextTextureCoords2 to TextureCoords1.
% The returned NewTextureCoords1 should generically have different dimensions in
% column, compared to that of TextureCoords1.
%

BBox = [-1,-1,1,1;-1,1,-1,1];
NextTextureCoords2 = [NextTextureCoords2,BBox];
TextureCoords1 = [TextureCoords1,BBox];
Tri = delaunayTriangulation(NextTextureCoords2');

NewTextureCoords1 = zeros(size(NextTextureCoords1));
nV = size(NewTextureCoords1,2);
nF = size(Tri.ConnectivityList,1);
for j = 1:nV
    if strcmpi(ProgressBar,'on')
        progressbar(j,nV,20);
    end
    BC = Tri.cartesianToBarycentric((1:nF)',repmat(NextTextureCoords1(:,j)',nF,1));
    tind = find(all(BC>-1e-10,2));
    if(numel(tind)>1)
        tind = tind(1);
    elseif numel(tind)==0
        warning(['Point ', num2str(j), ' was not found in any triangle. It stays where it was.']);
        NewTextureCoords1(:,j) = NextTextureCoords1(:,j);
        continue;
    end
    BC = BC(tind,:);
    NewTextureCoords1(:,j) = TextureCoords1(:,Tri.ConnectivityList(tind,:))*BC';
end

end

