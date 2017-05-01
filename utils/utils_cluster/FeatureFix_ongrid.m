function FeatureFix_ongrid(G1,G2,rslt_mat,TAXAind1,TAXAind2,LandmarksPath,MeshesPath,MeshSuffix,NumLandmark,options)

GM = load(G1);
GM = GM.G;
GN = load(G2);
GN = GN.G;

load(rslt_mat);
options.NumLandmark = str2double(NumLandmark);

tic;
disp(['Fixing features for ' GM.Aux.name ' vs ' GN.Aux.name '...']);
rslt = FeatureFix(GM,GN,TAXAind1,TAXAind2,options);
lk2 = GN.V(:,GetLandmarks(GN.Aux.name,LandmarksPath,[MeshesPath GN.Aux.name MeshSuffix],options));
lk1 = GN.V(:,rslt.ImprMap(GetLandmarks(GM.Aux.name,LandmarksPath,[MeshesPath GM.Aux.name MeshSuffix],options)));
% lk2 = GN.V(:,GetLandmarks(GN,LandmarksPath));
% lk1 = GN.V(:,rslt.ImprMap(GetLandmarks(GM,LandmarksPath)));
rslt.lkMSE = mean(sqrt(sum((lk2-lk1).^2)));
Imprrslt{TAXAind1,TAXAind2} = rslt;
save(rslt_mat,'Imprrslt');
disp(['Feature Fixing for ' GM.Aux.name ' vs ' GN.Aux.name ' done.']);
toc;

end

