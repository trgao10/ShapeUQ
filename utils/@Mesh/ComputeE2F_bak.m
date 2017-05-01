function E2F = ComputeE2F(G)

% ComputeE2F - compute faces adjacent to each edge
%
%   e2f = compute_edge_face_ring(faces);
%
%   e2f(i,j) and e2f(j,i) are the number of the two faces adjacent to
%   edge (i,j).
%
%   Copyright (c) 2007 Gabriel Peyre


% if (isfield(G, 'E2F') && ~isempty(G.E2F))
%     E2F = G.E2F;
%     return
% end

faces=G.F';

n = max(faces(:));
m = size(faces,2);
i = [faces(1,:) faces(2,:) faces(3,:)];
j = [faces(2,:) faces(3,:) faces(1,:)];
s = [1:m 1:m 1:m];

% first without duplicate
[~,I] = unique( i+(max(i)+1)*j );
% remaining items
J = setdiff(1:length(s), I);

% flip the duplicates
i1 = [i(I) j(J)];
j1 = [j(I) i(J)];
s = [s(I) s(J)];

% remove doublons
[~,I] = unique( i1+(max(i1)+1)*j1 );
i1 = i1(I); j1 = j1(I); s = s(I);

E2F = sparse(i1,j1,s,n,n);


% add missing points
I = find( E2F'~=0 );
I = I( E2F(I)==0 ); 
E2F( I ) = -1;
