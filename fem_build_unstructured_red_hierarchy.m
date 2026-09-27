function hierarchy = fem_build_unstructured_red_hierarchy(nBoundaryEdges, nInterior, nRefine, seed)

arguments
    nBoundaryEdges (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(nBoundaryEdges,2)} = 4
    nInterior (1,1) double {mustBeInteger,mustBePositive} = 10
    nRefine (1,1) double {mustBeInteger,mustBeNonnegative} = 4
    seed (1,1) double {mustBeInteger} = 20260821
end

[p,t,boundary] = coarse_mesh(nBoundaryEdges,nInterior,seed);
meshes = cell(nRefine+1,1);
P = cell(nRefine,1);
Pfull = cell(nRefine,1);
meshes{1} = finish_mesh(p,t,boundary);

for ell = 1:nRefine
    [p,t,boundary,Pfull{ell}] = red_refine( ...
        meshes{ell}.p,meshes{ell}.t,meshes{ell}.isBoundary);
    meshes{ell+1} = finish_mesh(p,t,boundary);
    P{ell} = Pfull{ell}(meshes{ell+1}.free,meshes{ell}.free);
end

hierarchy = struct('mesh',{meshes},'P',{P},'Pfull',{Pfull}, ...
    'seed',seed,'nBoundaryEdges',nBoundaryEdges,'nInterior',nInterior, ...
    'generator','constrained Delaunay coarse mesh + conforming red refinement');
end

function [p,t,boundary] = coarse_mesh(nb,nInterior,seed)
s = linspace(0,1,nb+1).';
bottom = [s(1:end-1),zeros(nb,1)];
right  = [ones(nb,1),s(1:end-1)];
top    = [flipud(s(2:end)),ones(nb,1)];
left   = [zeros(nb,1),flipud(s(2:end))];
bnd = [bottom;right;top;left];
nBoundary = size(bnd,1);
constraints = [(1:nBoundary).',[2:nBoundary,1].'];

oldState = rng;
cleanup = onCleanup(@() rng(oldState));
rng(seed,'twister');
interior = zeros(nInterior,2);
nAccepted = 0;
minSep = min(0.12,0.48/sqrt(nInterior));
for attempt = 1:200000
    if nAccepted == nInterior
        break
    end
    candidate = 0.07+0.86*rand(1,2);
    if nAccepted == 0 || all(vecnorm(interior(1:nAccepted,:)-candidate,2,2)>=minSep)
        nAccepted = nAccepted+1;
        interior(nAccepted,:) = candidate;
    end
end
if nAccepted ~= nInterior
    error('Could not generate the requested separated interior cloud.');
end

dt = delaunayTriangulation([bnd;interior],constraints);
p = dt.Points;
t = dt.ConnectivityList;
centroid = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;
inside = all(centroid>=-1e-13 & centroid<=1+1e-13,2);
t = orient_ccw(p,t(inside,:));
boundary = false(size(p,1),1);
boundary(1:nBoundary) = true;
clear cleanup
end

function mesh = finish_mesh(p,t,isBoundary)
t = orient_ccw(p,t);
isBoundary = logical(isBoundary(:));
free = find(~isBoundary);
boundary = find(isBoundary);
full2free = zeros(size(p,1),1);
full2free(free) = 1:numel(free);

p1 = p(t(:,1),:); p2 = p(t(:,2),:); p3 = p(t(:,3),:);
twiceArea = cross2(p2-p1,p3-p1);
edges = [vecnorm(p2-p1,2,2);vecnorm(p3-p2,2,2);vecnorm(p1-p3,2,2)];
a = vecnorm(p2-p3,2,2);
b = vecnorm(p1-p3,2,2);
c = vecnorm(p1-p2,2,2);
angles = [safe_acos((b.^2+c.^2-a.^2)./(2*b.*c)); ...
          safe_acos((a.^2+c.^2-b.^2)./(2*a.*c)); ...
          safe_acos((a.^2+b.^2-c.^2)./(2*a.*b))];

mesh = struct;
mesh.p = p;
mesh.t = t;
mesh.nodes = p;
mesh.elems = t;
mesh.free = free;
mesh.boundary = boundary;
mesh.isBoundary = isBoundary;
mesh.full2free = full2free;
mesh.hmax = max(edges);
mesh.hmin = min(edges);
mesh.minArea = min(twiceArea)/2;
mesh.minAngleDeg = min(angles)*180/pi;
mesh.nNode = size(p,1);
mesh.nElem = size(t,1);
mesh.nDof = numel(free);
end

function [pFine,tFine,boundaryFine,Pfull] = red_refine(p,t,boundary)
nOld = size(p,1);
nElem = size(t,1);
allEdges = [t(:,[1 2]);t(:,[2 3]);t(:,[3 1])];
[edges,~,edgeId] = unique(sort(allEdges,2),'rows');
nEdge = size(edges,1);
edgeCount = accumarray(edgeId,1,[nEdge,1]);
midpoint = (p(edges(:,1),:)+p(edges(:,2),:))/2;
pFine = [p;midpoint];
midId = nOld+(1:nEdge).';

e12 = edgeId(1:nElem);
e23 = edgeId(nElem+(1:nElem));
e31 = edgeId(2*nElem+(1:nElem));
m12 = midId(e12); m23 = midId(e23); m31 = midId(e31);
tFine = [t(:,1),m12,m31; ...
         m12,t(:,2),m23; ...
         m31,m23,t(:,3); ...
         m12,m23,m31];
tFine = orient_ccw(pFine,tFine);
boundaryFine = [boundary;edgeCount==1];

rows = [(1:nOld).';midId;midId];
cols = [(1:nOld).';edges(:,1);edges(:,2)];
vals = [ones(nOld,1);0.5*ones(2*nEdge,1)];
Pfull = sparse(rows,cols,vals,nOld+nEdge,nOld);
end

function t = orient_ccw(p,t)
p1 = p(t(:,1),:); p2 = p(t(:,2),:); p3 = p(t(:,3),:);
signedTwiceArea = cross2(p2-p1,p3-p1);
if any(abs(signedTwiceArea)<100*eps)
    error('Degenerate triangle detected.');
end
flip = signedTwiceArea<0;
t(flip,[2 3]) = t(flip,[3 2]);
end

function z = cross2(x,y)
z = x(:,1).*y(:,2)-x(:,2).*y(:,1);
end

function theta = safe_acos(x)
theta = acos(max(-1,min(1,x)));
end
