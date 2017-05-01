function [rslt] = PairlyImproveMaps(GM,GN,rslt12,rslt21,options)
%PAIRLYIMPROVEMAPS Summary of this function goes here
%   Detailed explanation goes here

Dist12 = rslt12.ImprDist;
Dist21 = rslt21.ImprDist;
VM12 = rslt12.TextureCoords1;
VN12 = rslt12.TextureCoords2;
VM21 = rslt21.TextureCoords2;
VN21 = rslt21.TextureCoords1;
ref12 = rslt12.ref;
ref21 = rslt21.ref;

while Dist12~=Dist21
    if Dist12<Dist21
        inv_map12 = knnsearch(VM12',VN12');
        [Dist21,~,VN21,VM21] = IterProjMoebius(GN,GM,inv_map12,ref21,options);
        if Dist12<Dist21
            break;
        end
    elseif Dist12>Dist21
        inv_map21 = knnsearch(VN21',VM21');
        [Dist12,~,VM12,VN12] = IterProjMoebius(GM,GN,inv_map21,ref12,options);
        if Dist12>Dist21
            break;
        end
    else
        break;
    end
end

options.GaussMinInds = 'on';
while Dist12~=Dist21
    if Dist12<Dist21
        inv_map12 = knnsearch(VM12',VN12');
        [Dist21,~,VN21,VM21] = IterProjMoebius(GN,GM,inv_map12,ref21,options);
        if Dist12<Dist21
            break;
        end
    elseif Dist12>Dist21
        inv_map21 = knnsearch(VN21',VM21');
        [Dist12,~,VM12,VN12] = IterProjMoebius(GM,GN,inv_map21,ref12,options);
        if Dist12>Dist21
            break;
        end
    else
        break;
    end
end

rslt.Gname1 = rslt12.Gname1;
rslt.Gname2 = rslt12.Gname2;
if Dist12<Dist21
    rslt.ImprDist = Dist12;
    rslt.TextureCoords1 = VM12;
    rslt.TextureCoords2 = VN12;
    rslt.ImprMap = knnsearch(VN12',VM12');
    rslt.invImprMap = knnsearch(VM12',VN12');
    rslt.ref = rslt12.ref;
else
    rslt.ImrpDist = Dist21;
    rslt.TextureCoords1 = VM21;
    rslt.TextureCoords2 = VN21;
    rslt.ImprMap = knnsearch(VN21',VM21');
    rslt.invImprMap = knnsearch(VM21',VN21');
    rslt.ref = rslt21.ref;
end


end

