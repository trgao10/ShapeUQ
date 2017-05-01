function TextureCoords = ResolveNaN(TextureCoords,OrigTextureCoords)

NaNInds = find(isnan(compl(TextureCoords)));
DelNaNTextureCoords = TextureCoords;
DelNaNTextureCoords(:,NaNInds) = [];
DelNaNTextureCoords(:,isnan(compl(DelNaNTextureCoords))) = [];
DelNaNTextureCoords(:,isinf(compl(DelNaNTextureCoords))) = [];
DelNaNTextureCoords_kdtree = kdtree_build(DelNaNTextureCoords');

if ~isempty(NaNInds)
    for j=1:length(NaNInds)
        progressbar(j,length(NaNInds),20);
        ThreeNN = kdtree_k_nearest_neighbors(DelNaNTextureCoords_kdtree,OrigTextureCoords(:,NaNInds(j)),3);
        TR = triangulation([1,2,3],OrigTextureCoords(:,ThreeNN)');
        BC = TR.cartesianToBarycentric(1,OrigTextureCoords(:,NaNInds(j))');
        TextureCoords(:,NaNInds(j)) = DelNaNTextureCoords(:,ThreeNN)*BC';
    end
end

end

