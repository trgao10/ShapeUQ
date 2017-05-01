function cPdist_ongrid(G1,G2,rslt_mat,TAXAind1,TAXAind2,LandmarksPath,MeshesPath,MeshSuffix,NumLandmark)

GM = load(G1);
GM = GM.G;
GN = load(G2);
GN = GN.G;

load(rslt_mat);

options.NumLandmark = str2double(NumLandmark);
options.FeatureType = 'ConfMax';
options.NumDensityPnts = 1000;
options.AngleIncrement = 0.01;
options.NumFeatureMatch = 4;
options.GaussMinMatch = 'off';
options.ProgressBar = 'off';

disp(['Comparing ' GM.Aux.name ' vs ' GN.Aux.name '...']);

rslt = GM.ComputeContinuousProcrustes(GN,options);
lk2 = GN.V(:,GetLandmarks(GN.Aux.name,LandmarksPath,[MeshesPath GN.Aux.name MeshSuffix],options));
lk1 = GN.V(:,rslt.cPmap(GetLandmarks(GM.Aux.name,LandmarksPath,[MeshesPath GM.Aux.name MeshSuffix],options)));
rslt.lkMSE = mean(sqrt(sum((lk2-lk1).^2)));

cPrslt{str2double(TAXAind1),str2double(TAXAind2)} = rslt;
save(rslt_mat,'cPrslt');

disp([GM.Aux.name ' vs ' GN.Aux.name ' done.']);

end

