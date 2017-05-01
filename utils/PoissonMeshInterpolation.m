function [reconMesh,ORT,PSD] = PoissonMeshInterpolation(domainMesh,MeshList,Weights,ORT,PSD)
%POISSONMESHINTERPOLATION Summary of this function goes here
%   Detailed explanation goes here

GroupSize = length(MeshList);
AffineTransformations = zeros(domainMesh.nF, GroupSize, 3, 2);

if (nargin<4)
    disp('Polar Decomposition...');
    ORT = zeros(domainMesh.nF, GroupSize, 3, 3);
    PSD = zeros(domainMesh.nF, GroupSize, 3, 3);
    
    cback = 0;
    for j=1:domainMesh.nF
        localF = domainMesh.F(:,j);
        localV = domainMesh.V(:,localF);
        DomainLocalFrame = [localV(1:2,2)-localV(1:2,1), localV(1:2,3)-localV(1:2,1)];
        for k=1:GroupSize
            localFk = MeshList{k}.F(:,j);
            localVk = MeshList{k}.V(:,localFk);
            kLocalFrame = [localVk(:,2)-localVk(:,1), localVk(:,3)-localVk(:,1)];
            
            AffineTransformations(j,k,:,:) = kLocalFrame/DomainLocalFrame;
            [U,S,V] = svd(squeeze(AffineTransformations(j,k,:,:)));
            S = [S,zeros(3,1)];
            makeRot = det(U)*det(V);
            if (makeRot > 0)
                V = [[V;zeros(1,2)],[0,0,1]'];
            else
                V = [[V;zeros(1,2)],[0,0,-1]'];
            end
            ORT(j,k,:,:) = U*V';
            PSD(j,k,:,:) = V*S*V';
        end
        
        for cc=1:cback
            fprintf('\b');
        end
        cback = fprintf(['%4d/' num2str(domainMesh.nF) ' done.\n'], j);
    end
end

%% reconstruction
%%%% interpolate affine transformation using polar decomposition
Weights = Weights/sum(Weights);

InterpAffineTransformations = zeros(domainMesh.nF, 3, 3);

cback = 0;
for j=1:domainMesh.nF
    InterpPSD = zeros(3,3);
%     InterpQua = zeros(4,1);
    if (GroupSize==2)
        InterpPSD = Weights(1)*squeeze(PSD(j,1,:,:))+Weights(2)*squeeze(PSD(j,2,:,:));
%         InterpQua = Weights(1)*qGetQ(squeeze(ORT(j,1,:,:)))+Weights(2)*qGetQ(squeeze(ORT(j,2,:,:)));
        InterpQua = slerp(qGetQ(squeeze(ORT(j,1,:,:))),qGetQ(squeeze(ORT(j,2,:,:))),1-Weights(1));
    else
        preInterpQua = zeros(GroupSize,4);
        for k=1:GroupSize
            InterpPSD = InterpPSD + Weights(k)*squeeze(PSD(j,k,:,:));
            preInterpQua(k,:) = qGetQ(squeeze(ORT(j,k,:,:)))';
%             InterpQua = InterpQua + Weights(k)*qGetQ(squeeze(ORT(j,k,:,:)));
        end
        InterpQua = wavg_quaternion_markley(preInterpQua, Weights');
    end
    InterpAffineTransformations(j, :, :) = qGetR(InterpQua/norm(InterpQua))*InterpPSD;
    
    for cc=1:cback
        fprintf('\b');
    end
    cback = fprintf(['%4d/' num2str(domainMesh.nF) ' done.\n'], j);
end

%%%% Poisson reconstruction
%%% compute divergence of the interpolated jacobian
[DivX,DivY] = domainMesh.ComputeDivergenceMatrix;
DivRhoX = DivX*squeeze(InterpAffineTransformations(:,1,1))...
    +DivY*squeeze(InterpAffineTransformations(:,1,2));
DivRhoY = DivX*squeeze(InterpAffineTransformations(:,2,1))...
    +DivY*squeeze(InterpAffineTransformations(:,2,2));
DivRhoZ = DivX*squeeze(InterpAffineTransformations(:,3,1))...
    +DivY*squeeze(InterpAffineTransformations(:,3,2));

%%% reconstruct
L = domainMesh.ComputeCotanLaplacian;
reconX = [L;ones(1,domainMesh.nV)]\[DivRhoX;0];
reconY = [L;ones(1,domainMesh.nV)]\[DivRhoY;0];
reconZ = [L;ones(1,domainMesh.nV)]\[DivRhoZ;0];

%%% check reconstructed mesh
reconMesh = Mesh('VF',[reconX';reconY';-reconZ'],domainMesh.F);

end

