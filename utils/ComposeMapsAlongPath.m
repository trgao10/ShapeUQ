function map12 = ComposeMapsAlongPath(InputPath,mapMatrix)
%COMPOSEMAPSALONGPATH Summary of this function goes here
%   Detailed explanation goes here

if length(InputPath)==1
    map12 = mapMatrix{InputPath(1), InputPath(1)};
    return;
end

map12 = mapMatrix{InputPath(1), InputPath(2)};

if (length(InputPath)>2)
    for j=3:length(InputPath)
        map12 = mapMatrix{InputPath(j-1), InputPath(j)}(map12);
    end
end

end

