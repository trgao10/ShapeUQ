function flatten_ongrid(mesh_file, sample_file)
%sample_relax_cP:   sample on one surface
%   Detailed explanation goes here

%==========================================================================
% Preprocessing
%==========================================================================

G = Mesh('off', mesh_file);
revName = strtok(mesh_file(end:-1:1),'/');
G.Aux.name = strtok(revName(end:-1:1),'_.');
[G.Aux.Area,G.Aux.Center] = G.Centralize('ScaleArea');
options.GaussMaxLocalWidth = 12; %% for Clement data set
options.GaussMinLocalWidth = 7; %% for Clement data set
G.ComputeMidEdgeUniformization(options); %%% default options only for PNAS

G.Nf = G.ComputeFaceNormals;
G.Nv = G.F2V'*G.Nf';
G.Nv = G.Nv'*diag(1./sqrt(sum((G.Nv').^2,1)));

%%% Compute cotangent Laplacian operator.
G.Aux.LB = G.ComputeCotanLaplacian;

%%% Save results to a .mat file.
save(sample_file, 'G');

end

