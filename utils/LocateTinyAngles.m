function [VertInds] = LocateTinyAngles(G,threshold)
%LOCATETINYANGLES Summary of this function goes here
%   Detailed explanation goes here

Angles = G.ComputeTriangleAngles;

[r,c] = find(Angles<threshold);
VertInds = [];
for j=1:length(r)
    VertInds = [VertInds G.F(c(j),r(j))];
end

end

