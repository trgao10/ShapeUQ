function Imprdist_ongrid(G1,G2,rslt_mat,TAXAind1,TAXAind2,LandmarksPath,MeshesPath,MeshSuffix,NumLandmark,ImprType,FeatureFix,cPdistMatrix,cPmapsMatrix,options)

GM = load(G1);
GM = GM.G;
GN = load(G2);
GN = GN.G;

load(rslt_mat);

options.NumLandmark = str2double(NumLandmark);
options.ImprType = ImprType;
options.FeatureFix = FeatureFix;
options.Angle = 0.25; % ViterbiAngle
options.alpha = 'auto';    % actually this doesn't matter for
                           % cPComposedLAST; for setting it up
                           % on grid, it suffices if one
                           % specifie"cPLASTPath" in "options".
% options.alpha = 1+sqrt(2); % LAST/ComposedLAST; scalar>1 or 'auto'
options.SmoothMap = 1;
options.ProgressBar = 'off';

tic;
disp(['Improving ' GM.Aux.name ' vs ' GN.Aux.name ' using ' ImprType '...']);
rslt = GM.ImproveMap(GN,cPdistMatrix,cPmapsMatrix,options.TaxaCode,options);
lk2 = GN.V(:,GetLandmarks(GN.Aux.name,LandmarksPath,[MeshesPath GN.Aux.name MeshSuffix],options));
lk1 = GN.V(:,rslt.ImprMap(GetLandmarks(GM.Aux.name,LandmarksPath,[MeshesPath GM.Aux.name MeshSuffix],options)));
% lk2 = GN.V(:,GetLandmarks(GN,LandmarksPath));
% lk1 = GN.V(:,rslt.ImprMap(GetLandmarks(GM,LandmarksPath)));
rslt.lkMSE = mean(sqrt(sum((lk2-lk1).^2)));
Imprrslt{TAXAind1,TAXAind2} = rslt;
save(rslt_mat,'Imprrslt');
disp([GM.Aux.name ' vs ' GN.Aux.name ' done.']);
toc;

end

