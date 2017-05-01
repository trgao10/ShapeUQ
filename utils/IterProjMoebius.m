function [Dist12,proj_map12,VM12,VN12] = IterProjMoebius(GM,GN,map12,ref12,options)
%ITERPROJMOEBIUS: Iteratively Find the Moebius transform from s1 to s2 that
%                 is closest to map12. In some sense this is "projecting a
%                 map to the subspace of Moebius transforms
%   GM:           triangular mesh with M vertices
%   GN:           triangular mesh with N vertices
%   map12:        Mx1 vector, max(map12)<=N, so that GN.V(:,map12) becomes
%                 a vector of size Mx1, which can be used as the texture
%                 coordinates for GM
%   ref12:        0 (no reflection) or 1 (reflection)
%
%   Tingran Gao, Duke University
%   trgao10@math.duke.edu
%

[VM12,VN12,proj_map12] = ProjMoebius(GM,GN,map12,ref12,options);
Dist12 = MapToDist(GM.V,GN.V,proj_map12,GM.Aux.VertArea);
disp(['Vertex permutation with cPValue ' num2str(MapToDist(GM.V,GN.V,map12,GM.Aux.VertArea))...
      ' projected to a map cPValue ' num2str(Dist12)]);

[VM_iter,VN_iter,proj_map12_iter] = ProjMoebius(GM,GN,proj_map12,ref12,options);
while Dist12>MapToDist(GM.V,GN.V,proj_map12_iter,GM.Aux.VertArea)
    disp(['Decreased cPValue from ' num2str(Dist12) ' to ' num2str(MapToDist(GM.V,GN.V,proj_map12_iter,GM.Aux.VertArea))]);
    Dist12 = MapToDist(GM.V,GN.V,proj_map12_iter,GM.Aux.VertArea);
    proj_map12 = proj_map12_iter;
    VM12 = VM_iter;
    VN12 = VN_iter;
    [VM_iter,VN_iter,proj_map12_iter] = ProjMoebius(GM,GN,proj_map12,ref12,options);
end

options.GaussMinMatch = 'on';
[VM_iter,VN_iter,proj_map12_iter] = ProjMoebius(GM,GN,proj_map12,ref12,options);
while Dist12>MapToDist(GM.V,GN.V,proj_map12_iter,GM.Aux.VertArea)
    disp(['Matching GaussMinInds decreased cPValue from ' num2str(Dist12) ' to ' num2str(MapToDist(GM.V,GN.V,proj_map12_iter,GM.Aux.VertArea))]);
    Dist12 = MapToDist(GM.V,GN.V,proj_map12_iter,GM.Aux.VertArea);
    proj_map12 = proj_map12_iter;
    VM12 = VM_iter;
    VN12 = VN_iter;
    [VM_iter,VN_iter,proj_map12_iter] = ProjMoebius(GM,GN,proj_map12,ref12,options);
end

disp('Finished iterative Moebius projection.');

end

