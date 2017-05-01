function V=spin(T,V,rho)
% input: T (3xN) and V (3xM) so that
%  trimesh(T',V(1,:),V(2,:),V(3,:)) plots the triangular mesh,
%  rho (1xN) values of conformal scaling
% output: V (3xM) vertices of transformed mesh
%
% Author: Jan Hakenberg
% http://www.hakenberg.de/diffgeo/differential_geometry.htm#Spin
%


nT=size(T,2);
nV=size(V,2);
plc=-3:0;

E=sparse(4*nV,4*nV);
edg=zeros(4,3);
for c1=1:nT
  tri=T(:,c1);
  pnt=V(:,tri);
  A=norm( cross( pnt(:,2)-pnt(:,1) , pnt(:,3)-pnt(:,1) ) )/2;
  a=-1/(4*A);
  b=rho(c1)/6;
  c=jiH([A*rho(c1)*rho(c1)/9 0 0 0]);
  for c2=1:3
    edg(:,c2)=[0;V(:, tri(mod(c2+1,3)+1) )-V(:, tri(mod(c2+0,3)+1) )];
  end  
  ini=[tri(1)*4+plc tri(2)*4+plc tri(3)*4+plc];
  E(ini,ini)=E(ini,ini) + [
    jiH(a*jiH(edg(:,1))*edg(:,1)) + c ...
    jiH(a*jiH(edg(:,1))*edg(:,2) + b*(edg(:,2)-edg(:,1)))+c ...
    jiH(a*jiH(edg(:,1))*edg(:,3) + b*(edg(:,3)-edg(:,1)))+c
    jiH(a*jiH(edg(:,2))*edg(:,1) + b*(edg(:,1)-edg(:,2)))+c ...
    jiH(a*jiH(edg(:,2))*edg(:,2)) + c ...
    jiH(a*jiH(edg(:,2))*edg(:,3) + b*(edg(:,3)-edg(:,2)))+c
    jiH(a*jiH(edg(:,3))*edg(:,1) + b*(edg(:,1)-edg(:,3)))+c ...
    jiH(a*jiH(edg(:,3))*edg(:,2) + b*(edg(:,2)-edg(:,3)))+c ...
    jiH(a*jiH(edg(:,3))*edg(:,3)) + c ];
  if ~mod(c1,500); fprintf('.'); end
end
fprintf('\n')

lam=zeros(4*nV,1);
lam(1:4:end)=1;  
for c1=1:11
  cnv=lam;
  lam=E\lam;
  lam=lam/norm(lam);
end
res=(E*lam)./lam;  
fprintf('mean %e, var %e, delta %e\n',mean(res),var(res),norm(cnv-lam))

L  =sparse(4*nV,4*nV);
ome=zeros(4*nV,1);
for c1=1:nT
  for c2=1:3
    k0=T(mod(c2-1,3)+1,c1);
    k1=T(mod(c2+0,3)+1,c1);
    k2=T(mod(c2+1,3)+1,c1);
    u1=V(:,k1)-V(:,k0);
    u2=V(:,k2)-V(:,k0);
    cta=dot(u1,u2) / norm( cross(u1,u2) );
    h=jiH([cta*0.5 0 0 0]);
    ini=[k1*4+plc  k2*4+plc];
    L(ini,ini)=L(ini,ini)+[ h -h;-h h];
    if k1>k2
      k3=k1; k1=k2; k2=k3; % swap
    end
    lm1=jiH(lam(k1*4+plc));
    lm2=jiH(lam(k2*4+plc));
    edv=jiH([0;V(:,k2)-V(:,k1)]);
    til=lm1'*edv*lm1/3 + lm1'*edv*lm2/6 + lm2'*edv*lm1/6 + lm2'*edv*lm2/3;
    ome(k1*4+plc,1)=ome(k1*4+plc,1)-cta*til(:,1)/2;
    ome(k2*4+plc,1)=ome(k2*4+plc,1)+cta*til(:,1)/2;
  end
  if ~mod(c1,500); fprintf('.'); end
end
fprintf('\n')

ome=reshape(ome,[4 nV]);
ome=ome-repmat(mean(ome,2),[1 nV]);
ome=reshape(ome,[4*nV 1]);
ome=L\ome;
ome=reshape(ome,[4 nV]);
ome=ome-repmat(mean(ome,2),[1 nV]);
nrm=sum(ome.*ome,1);
ome=ome/sqrt(max(nrm));
V=ome(2:end,:);

function h=jiH(pnt)
a=pnt(1); b=pnt(2); c=pnt(3); d=pnt(4);
h=[ a -b -c -d
    b  a -d  c
    c  d  a -b
    d -c  b  a];
