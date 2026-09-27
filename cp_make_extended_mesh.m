function mesh = cp_make_extended_mesh(n, deltaMax, seed)

arguments
    n (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(n,2)}
    deltaMax (1,1) double {mustBeNonnegative}
    seed (1,1) double {mustBeInteger} = 20260821
end

h = 1/n;
nExt = max(1, ceil(deltaMax/h) + 1);
idx = (-nExt):(n+nExt);
xv = idx*h;
[X,Y] = meshgrid(xv,xv);
p = [X(:),Y(:)];

old = rng;
clean = onCleanup(@() rng(old));
rng(seed,'twister');
tol = 100*eps;
movable = p(:,1)>tol & p(:,1)<1-tol & p(:,2)>tol & p(:,2)<1-tol;
jitter = 0.16*h*(2*rand(nnz(movable),2)-1);
p(movable,:) = p(movable,:) + jitter;

nx = numel(xv);
t = zeros(2*(nx-1)^2,3);
region = false(2*(nx-1)^2,1);
e = 0;
for j = 1:nx-1
    for i = 1:nx-1
        bl = i + nx*(j-1);
        br = bl+1;
        tl = bl+nx;
        tr = tl+1;
        if rand < 0.5
            loc = [bl br tr; bl tr tl];
        else
            loc = [bl br tl; br tr tl];
        end
        for k = 1:2
            e = e+1;
            t(e,:) = loc(k,:);
            c = mean(p(loc(k,:),:),1);
            region(e) = c(1)>0 && c(1)<1 && c(2)>0 && c(2)<1;
        end
    end
end
clear clean

t = orient_ccw(p,t);
onBoundary = abs(p(:,1))<tol | abs(p(:,1)-1)<tol | ...
             abs(p(:,2))<tol | abs(p(:,2)-1)<tol;
insideNode = p(:,1)>tol & p(:,1)<1-tol & p(:,2)>tol & p(:,2)<1-tol;
free = find(insideNode);
full2free = zeros(size(p,1),1);
full2free(free) = 1:numel(free);

triP = reshape(p(t',:),3,[],2);
centres = squeeze(mean(triP,1));
bbox = [min(reshape(triP(:,:,1),3,[]),[],1)', ...
        max(reshape(triP(:,:,1),3,[]),[],1)', ...
        min(reshape(triP(:,:,2),3,[]),[],1)', ...
        max(reshape(triP(:,:,2),3,[]),[],1)'];

mesh = struct;
mesh.p = p;
mesh.t = t;
mesh.nodes = p;
mesh.elems = t;
mesh.isInteriorElem = region;
mesh.free = free;
mesh.boundary = onBoundary;
mesh.full2free = full2free;
mesh.centres = centres;
mesh.bbox = bbox;
mesh.h = h;
mesh.hmax = max_triangle_edge(p,t(region,:));
mesh.deltaMax = deltaMax;
mesh.seed = seed;
mesh.omega = [0 1 0 1];
end

function t = orient_ccw(p,t)
a = p(t(:,1),:); b = p(t(:,2),:); c = p(t(:,3),:);
d = (b(:,1)-a(:,1)).*(c(:,2)-a(:,2)) - ...
    (b(:,2)-a(:,2)).*(c(:,1)-a(:,1));
if any(abs(d)<1e-14)
    error('cp_make_extended_mesh:degenerate','Degenerate triangle.');
end
flip = d<0;
t(flip,[2 3]) = t(flip,[3 2]);
end

function hmax = max_triangle_edge(p,t)
if isempty(t), hmax=NaN; return; end
d1 = vecnorm(p(t(:,1),:)-p(t(:,2),:),2,2);
d2 = vecnorm(p(t(:,2),:)-p(t(:,3),:),2,2);
d3 = vecnorm(p(t(:,3),:)-p(t(:,1),:),2,2);
hmax = max([d1;d2;d3]);
end
