function pc = cp_paircell_4d(T,S,Z,delta,opts)
arguments
    T (3,2) double
    S (3,2) double
    Z (3,2) double
    delta (1,1) double {mustBePositive}
    opts.tol (1,1) double {mustBePositive} = 5e-11
end
[AT,bT]=cp_triangle_halfspace(T);
[AS,bS]=cp_triangle_halfspace(S);
[AZ,bZ]=cp_triangle_halfspace(Z);
A=[AT,zeros(3,2); zeros(3,2),AS; -AZ/delta,AZ/delta];
b=[bT;bS;bZ];

comb=nchoosek(1:size(A,1),4);
V=zeros(size(comb,1),4); nv=0;
scale=max(1,max(abs([T(:);S(:);delta*Z(:)])));
tol=opts.tol*scale;
for k=1:size(comb,1)
    Ak=A(comb(k,:),:);
    if rcond(Ak)<1e-12, continue; end
    x=Ak\b(comb(k,:));
    if all(A*x<=b+tol)
        nv=nv+1; V(nv,:)=x';
    end
end
V=V(1:nv,:);
if nv<5
    pc=empty_pc(V,A,b); return
end
V=uniquetol(V,max(1e-12,20*opts.tol),'ByRows',true,'DataScale',scale);
if size(V,1)<5 || rank(V(2:end,:)-V(1,:),1e-11*scale)<4
    pc=empty_pc(V,A,b); return
end

if size(V,1)==5
    simp=1:5;
else
    try
        simp=delaunayn(V,{'Qbb','Qc','Qz','Qt'});
    catch
        simp=delaunayn(V,{'QJ','Qbb','Qc'});
    end
end

ns=size(simp,1);
vol=zeros(ns,1); keep=false(ns,1);
for k=1:ns
    R=V(simp(k,:),:);
    vol(k)=abs(det(R(2:5,:)-R(1,:)))/24;
    keep(k)=vol(k)>max(1e-15,opts.tol^2)*scale^4;
end
simp=simp(keep,:); vol=vol(keep);
if isempty(simp), pc=empty_pc(V,A,b); return; end

n=4;
aa=(1+n/sqrt(n+2))/(n+1);
bb=(1-1/sqrt(n+2))/(n+1);
nq=5*size(simp,1);
q=zeros(nq,4); w=zeros(nq,1); owner=zeros(nq,1);
pos=0;
for k=1:size(simp,1)
    R=V(simp(k,:),:);
    total=sum(R,1);
    for i=1:5
        pos=pos+1;
        q(pos,:)=aa*R(i,:)+bb*(total-R(i,:));
        w(pos)=vol(k)/5;
        owner(pos)=k;
    end
end
pc=struct('vertices',V,'simplices',simp,'simplexVolume',vol, ...
    'volume',sum(vol),'q',q,'w',w,'simplexOwner',owner,'A',A,'b',b, ...
    'nVertices',size(V,1),'nSimplices',size(simp,1),'empty',false);
end

function pc=empty_pc(V,A,b)
pc=struct('vertices',V,'simplices',zeros(0,5),'simplexVolume',zeros(0,1), ...
    'volume',0,'q',zeros(0,4),'w',zeros(0,1),'simplexOwner',zeros(0,1), ...
    'A',A,'b',b,'nVertices',size(V,1),'nSimplices',0,'empty',true);
end
