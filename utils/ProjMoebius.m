function [V1,V2,proj_map12] = ProjMoebius(GM,GN,map12,ref12,options)
%PROJMOEBIUS: Find the Moebius transform from s1 to s2 that is closest to
%             map12. In some sense this is "projecting a map to the
%             subspace of Moebius transforms
%   GM:       triangular mesh with M vertices
%   GN:       triangular mesh with N vertices
%   map12:    Mx1 vector, max(map12)<=N, so that GN.V(:,map12) becomes a
%             vector of size Mx1, which can be used as the texture
%             coordinates for GM
%   ref12:    0 (no reflection) or 1 (reflection)
%
%   Tingran Gao, Duke University
%   trgao10@math.duke.edu
%

if nargin<5
    options = [];
end
GaussMinMatch = getoptions(options,'GaussMinMatch','off');

%%% check for NaN's in the uniformization of GM
ts = GM.Aux.UniformizationV(1,:)+1i*GM.Aux.UniformizationV(2,:);
delInds = isnan(ts);
ts(delInds) = [];
%%% check for NaN's in the uniformization of GM
V2 = GN.Aux.UniformizationV(1:2,:);
V2(:,isnan(compl(V2))) = ones(2,sum(isnan(compl(V2))));
ref_GM = V2(1:2,map12);
ref_GM(:,delInds) = [];

if ref12==1
    V2(2,:) = -V2(2,:);
    ref_GM(2,:) = -ref_GM(2,:);
end

VertArea = GM.Aux.VertArea';

FeaturesM = GM.Aux.ConfMaxInds;
FeaturesN = GN.Aux.ConfMaxInds;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% starts Moebius projection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
best_err = Inf;
best_a = 0;
best_tet = 0;

for jj=1:length(FeaturesM)
    progressbar(jj,length(FeaturesM),10);
    for kk=1:length(FeaturesN)
        z_0 = GM.Aux.UniformizationV(1,FeaturesM(jj))+1i*GM.Aux.UniformizationV(2,FeaturesM(jj));
        w_0 = V2(1,FeaturesN(kk))+1i*V2(2,FeaturesN(kk));
        
        for tet=0:0.01:2*pi %traverse angles: sould be in range 0.05-0.1
            [a] = CORR_evaluate_disc_moebius_from_tet(tet,z_0,w_0);
            if(a*conj(a) > 0.9999)
                err = Inf;
            else
                % Push GM to GN by m
                m = [exp(1i*tet) -a*exp(1i*tet); -conj(a) 1];%takes z_0 -> w_0
                push_GM = CORR_apply_moebius_as_matrix(m,ts);
                push_GM = [real(push_GM);imag(push_GM)];
                err = sum((push_GM-ref_GM).^2)*VertArea;
            end
            % Record if best so far
            if (err < best_err)
                best_err = err;
                best_a = a;
                best_tet = tet;
            end
        end
    end
end

m = [exp(1i*best_tet) -best_a*exp(1i*best_tet); -conj(best_a) 1];
pushGM = CORR_apply_moebius_as_matrix(m,compl(GM.Aux.UniformizationV));
pushGM(isnan(pushGM)) = 1+1i;
V2 = GN.Aux.UniformizationV(1:2,:);
V2(:,isnan(compl(V2))) = ones(2,sum(isnan(compl(V2))));
if ref12==1
    V2(2,:) = -V2(2,:);
end
TextureCoords2_kdtree = kdtree_build(V2');

%%% match features in Euclidean space
[~,TPS_EUC_FEATURESN,preTPS_EUC_FEATURESM] = FindEuclideanMutuallyNearestNeighbors(GM,GN,map12,'ConfMax');

TPS_EUC_FEATURE_MATCH_M = DISCtoPLANE([real(pushGM(preTPS_EUC_FEATURESM));imag(pushGM(preTPS_EUC_FEATURESM))]','d2p');
TPS_EUC_FEATURE_MATCH_N = DISCtoPLANE(V2(:,TPS_EUC_FEATURESN)','d2p');

TPS_FEATURESM = TPS_EUC_FEATURE_MATCH_M;
TPS_FEATURESN = TPS_EUC_FEATURE_MATCH_N;

if (length(TPS_FEATURESM)>3) % TPS (Thin Plate Spline)
    tP = DISCtoPLANE([real(pushGM);imag(pushGM)]','d2p');
    [ftps] = TEETH_calc_tps(TPS_FEATURESM,TPS_FEATURESN-TPS_FEATURESM);
    pt = tP + TEETH_eval_tps(ftps,tP);
    V1 = DISCtoPLANE(pt,'p2d')';
elseif (length(TPS_FEATURESM)==3) % Affine Transformation
    tP = DISCtoPLANE([real(pushGM);imag(pushGM)]','d2p');
    [A,b] = PlanarThreePtsDeform(TPS_FEATURESM,TPS_FEATURESN);
    pt = [A,b]*[tP';ones(1,size(tP,1))];
    V1 = DISCtoPLANE(pt','p2d')';
end
proj_map12 = kdtree_nearest_neighbor(TextureCoords2_kdtree, V1');

if strcmpi(GaussMinMatch,'on')
    disp('Matching GaussMinInds');
    [~,InterpGaussMinInds2,preInterpGaussMinInds1] = FindEuclideanMutuallyNearestNeighbors(GM,GN,proj_map12,'GaussMin');
    TPS_GaussMinCoords1 = DISCtoPLANE([real(pushGM(preInterpGaussMinInds1));imag(pushGM(preInterpGaussMinInds1))]','d2p');
    TPS_GaussMinCoords2 = DISCtoPLANE(V2(:,InterpGaussMinInds2)','d2p');
    TPS_FEATURESM = [TPS_FEATURESM;TPS_GaussMinCoords1];
    TPS_FEATURESN = [TPS_FEATURESN;TPS_GaussMinCoords2];
    if (length(TPS_FEATURESM)>3) % TPS (Thin Plate Spline)
        tP = DISCtoPLANE([real(pushGM);imag(pushGM)]','d2p');
        [ftps] = TEETH_calc_tps(TPS_FEATURESM,TPS_FEATURESN-TPS_FEATURESM);
        pt = tP + TEETH_eval_tps(ftps,tP);
        V1 = DISCtoPLANE(pt,'p2d')';
    elseif (length(TPS_FEATURESM)==3) % Affine Transformation
        tP = DISCtoPLANE([real(pushGM);imag(pushGM)]','d2p');
        [A,b] = PlanarThreePtsDeform(TPS_FEATURESM,TPS_FEATURESN);
        pt = [A,b]*[tP';ones(1,size(tP,1))];
        V1 = DISCtoPLANE(pt','p2d')';
    end
    proj_map12 = kdtree_nearest_neighbor(TextureCoords2_kdtree, V1');
    disp('Done.');
end

if ref12==1
    V1(2,:) = -V1(2,:);
    V2(2,:) = -V2(2,:);
end

end
