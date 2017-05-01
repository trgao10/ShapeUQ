function [min_path,map12] = FindViterbiPath(GM,GN,DistMatrix,MapMatrix,TaxaCode,options)
%VIEWRMALONGVITERBI Summary of this function goes here
%   Detailed explanation goes here

if nargin<6
    options = [];
end
if ~isfield(options,'R')
    warning('Should include alignment R in "options"');
end
if ~isfield(options,'T')
    warning('Should include path length upper bound T in "options"');
end

R = getoptions(options,'R',eye(3));
T = getoptions(options,'T',30);
Angle = getoptions(options,'Angle','off');
if ~strcmpi(Angle,'off')
    disp('Angle Costs Added.');
end
tGMV = R*GM.V;
GroupSize = length(DistMatrix);
TAXAind1 = find(strcmpi(TaxaCode,GM.Aux.name));
TAXAind2 = find(strcmpi(TaxaCode,GN.Aux.name));

if TAXAind1==TAXAind2
    min_path = TAXAind1;
    map12 = MapMatrix{TAXAind1,TAXAind2};
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IntermediatePaths = cell(1,GroupSize);
IntermediatePathProbs = zeros(1,GroupSize);
IntermediateTempProbs = zeros(1,GroupSize);
FinalPaths = cell(1,T);

ReverseIntermediatePaths = cell(1,GroupSize);
ReverseIntermediatePathProbs = zeros(1,GroupSize);
ReverseIntermediateTempProbs = zeros(1,GroupSize);
ReverseFinalPaths = cell(1,T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recursively find most likely paths with number of edges from 1 to T
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['Find paths with number of edges from 1 to ' num2str(T) ' ...']);
for j=1:T
    progressbar(j,T);
    if j==1 % initialize data structures
        for k=1:GroupSize
            IntermediatePaths{k} = [TAXAind1,k];
            ReverseIntermediatePaths{k} = [TAXAind2,k];
        end
        TransProbs = DistMatrix(TAXAind1,:);
        TransProbs(TAXAind1) = max(TransProbs);
        TransProbs = max(TransProbs)-TransProbs;
        IntermediatePathProbs = TransProbs/sum(TransProbs);
        FinalPaths{j} = IntermediatePaths{TAXAind2};
        
        TransProbs = DistMatrix(TAXAind2,:);
        TransProbs(TAXAind2) = max(TransProbs);
        TransProbs = max(TransProbs)-TransProbs;
        ReverseIntermediatePathProbs = TransProbs/sum(TransProbs);
        ReverseFinalPaths{j} = ReverseIntermediatePaths{TAXAind1};
    else % update according to history
        NewIntermediatePaths = cell(size(IntermediatePaths));
        NewIntermediatePathProbs = zeros(size(IntermediatePathProbs));
        ReverseNewIntermediatePaths = cell(size(ReverseIntermediatePaths));
        ReverseNewIntermediatePathProbs = zeros(size(ReverseIntermediatePathProbs));
        for k=1:GroupSize
            % connect each intermediate path to k, find the best path that
            % connects to k, update IntermediatePaths{k}
            for kk=1:GroupSize
                % CurrentPath connects TAXAind1 to kk
                CurrentPath = IntermediatePaths{kk};
                % CurrentReversePath connects TAXAind2 to kk
                ReverseCurrentPath = ReverseIntermediatePaths{kk};
                
                % DistCosts: from kk to all other nodes in one step
                DistCosts = DistMatrix(kk,:);
                DistCosts(CurrentPath) = max(DistCosts); % avoid cycles
                DistCosts = max(DistCosts)-DistCosts;
                DistCosts = DistCosts/sum(DistCosts);
                
                ReverseDistCosts = DistMatrix(kk,:);
                ReverseDistCosts(ReverseCurrentPath) = max(ReverseDistCosts); % avoid cycles
                ReverseDistCosts = max(ReverseDistCosts)-ReverseDistCosts;
                ReverseDistCosts = ReverseDistCosts/sum(ReverseDistCosts);
                
                if ~strcmpi(Angle,'off')
                    % AngleCosts: from kk to all other nodes in one step
                    % similar for ReverseAngleCosts
                    AngleCosts = zeros(size(DistCosts));
                    ReverseAngleCosts = zeros(size(ReverseDistCosts));
                    for tt=1:GroupSize
                        Node = kk;
                        TipJ = CurrentPath(end-1);
                        ReverseTipJ = ReverseCurrentPath(end-1);
                        TipK = tt;
                        if (Node==TipJ)||(Node==TipK)||(TipJ==TipK)
                            AngleCosts(tt) = 1;
                        else
                            DistJK = DistMatrix(TipJ,TipK);
                            DistNJ = DistMatrix(Node,TipJ);
                            DistNK = DistMatrix(Node,TipK);
                            AngleCosts(tt) = (DistNJ^2+DistNK^2-DistJK^2)/(2*DistNJ*DistNK);
                        end
                        if (Node==ReverseTipJ)||(Node==TipK)||(ReverseTipJ==TipK)
                            ReverseAngleCosts(tt) = 1;
                        else
                            DistJK = DistMatrix(ReverseTipJ,TipK);
                            DistNJ = DistMatrix(Node,ReverseTipJ);
                            DistNK = DistMatrix(Node,TipK);
                            ReverseAngleCosts(tt) = (DistNJ^2+DistNK^2-DistJK^2)/(2*DistNJ*DistNK);
                        end
                    end
                    ReverseAngleCosts(ReverseCurrentPath) = 1; % avoid cycles
                    ReverseAngleCosts = 1-ReverseAngleCosts;
                    ReverseAngleCosts = ReverseAngleCosts/sum(ReverseAngleCosts);
                    
                    TransProbs = Angle*DistCosts+(1-Angle)*AngleCosts;
                    TransProbs = TransProbs/sum(TransProbs);
                    IntermediateTempProbs(kk) = IntermediatePathProbs(kk)*TransProbs(k);
                    
                    TransProbs = Angle*ReverseDistCosts+(1-Angle)*ReverseAngleCosts;
                    TransProbs = TransProbs/sum(TransProbs);
                    ReverseIntermediateTempProbs(kk) = ReverseIntermediatePathProbs(kk)*TransProbs(k);
                else
                    TransProbs = DistCosts;
                    TransProbs = TransProbs/sum(TransProbs);
                    IntermediateTempProbs(kk) = IntermediatePathProbs(kk)*TransProbs(k);
                    
                    TransProbs = ReverseDistCosts;
                    TransProbs = TransProbs/sum(TransProbs);
                    ReverseIntermediateTempProbs(kk) = ReverseIntermediatePathProbs(kk)*TransProbs(k);
                end
            end
            [NewIntermediatePathProbs(k),ind] = max(IntermediateTempProbs);
            NewIntermediatePaths{k} = [IntermediatePaths{ind},k];
            
            [ReverseNewIntermediatePathProbs(k),ind] = max(ReverseIntermediateTempProbs);
            ReverseNewIntermediatePaths{k} = [ReverseIntermediatePaths{ind},k];
        end
        IntermediatePaths = NewIntermediatePaths;
        IntermediatePathProbs = NewIntermediatePathProbs;
        FinalPaths{j} = IntermediatePaths{TAXAind2};
        
        ReverseIntermediatePaths = ReverseNewIntermediatePaths;
        ReverseIntermediatePathProbs = ReverseNewIntermediatePathProbs;
        ReverseFinalPaths{j} = ReverseIntermediatePaths{TAXAind1};
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the path from FinalPaths with lowest CP functional value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Evaluate maps along all optimal paths ...');
CandidateCPValues = zeros(1,T);
ReverseCandidateCPValues = zeros(1,T);
for j=1:T
    progressbar(j,T);
    map12 = ComposeMapsAlongPath(FinalPaths{j},MapMatrix);
    map21 = ComposeMapsAlongPath(ReverseFinalPaths{j},MapMatrix);
    CandidateCPValues(j) = MapToDist(tGMV,GN.V,map12,GM.Aux.VertArea);
    ReverseCandidateCPValues(j) = MapToDist(GN.V,tGMV,map21,GN.Aux.VertArea);
end

[~,ind] = min(CandidateCPValues);
ViterbiPath = FinalPaths{ind};
[~,ind] = min(ReverseCandidateCPValues);
ReverseViterbiPath = ReverseFinalPaths{ind};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare four different maps composed along each of the fourViterbi paths
% pick the best one
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ViterbiMap12 = ComposeMapsAlongPath(ViterbiPath,MapMatrix);
ViterbiMap12CPValue = MapToDist(tGMV,GN.V,ViterbiMap12,GM.Aux.VertArea);
ViterbiMap21 = ComposeMapsAlongPath(ViterbiPath(end:-1:1),MapMatrix);
ViterbiMap21CPValue = MapToDist(GN.V,tGMV,ViterbiMap21,GN.Aux.VertArea);
ReverseViterbiMap12 = ComposeMapsAlongPath(ReverseViterbiPath(end:-1:1),MapMatrix);
ReverseViterbiMap12CPValue = MapToDist(tGMV,GN.V,ReverseViterbiMap12,GM.Aux.VertArea);
ReverseViterbiMap21 = ComposeMapsAlongPath(ReverseViterbiPath,MapMatrix);
ReverseViterbiMap21CPValue = MapToDist(GN.V,tGMV,ReverseViterbiMap21,GN.Aux.VertArea);

[~,MinInd] = min([ViterbiMap12CPValue,ViterbiMap21CPValue,ReverseViterbiMap12CPValue,ReverseViterbiMap21CPValue]);
switch MinInd
    case 1
        map12 = ViterbiMap12;
        map21 = knnsearch(GN.Aux.UniformizationV(1:2,map12)',GN.Aux.UniformizationV(1:2,:)');
        min_path = ViterbiPath;
    case 2
        map21 = ViterbiMap21;
        map12 = knnsearch(GM.Aux.UniformizationV(1:2,map21)',GM.Aux.UniformizationV(1:2,:)');
        min_path = ViterbiPath;
    case 3
        map12 = ReverseViterbiMap12;
        map21 = knnsearch(GN.Aux.UniformizationV(1:2,map12)',GN.Aux.UniformizationV(1:2,:)');
        min_path = ReverseViterbiPath(end:-1:1);
    case 4
        map21 = ReverseViterbiMap21;
        map12 = knnsearch(GM.Aux.UniformizationV(1:2,map21)',GM.Aux.UniformizationV(1:2,:)');
        min_path = ReverseViterbiPath(end:-1:1);
end

disp('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
disp([GM.Aux.name ' to ' GN.Aux.name]);
disp('Viterbi Path:');
disp(TaxaCode(min_path));
disp(['Viterbi CPValue: ' num2str(MapToDist(tGMV,GN.V,map12,GM.Aux.VertArea))]);

end

