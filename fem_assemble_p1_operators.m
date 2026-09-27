function fem = fem_assemble_p1_operators(mesh)

if isfield(mesh,'p'), p = mesh.p; else, p = mesh.nodes; end
if isfield(mesh,'t'), t = mesh.t; else, t = mesh.elems; end
free = mesh.free(:);
nNode = size(p,1);
nElem = size(t,1);

x1 = p(t(:,1),1); y1 = p(t(:,1),2);
x2 = p(t(:,2),1); y2 = p(t(:,2),2);
x3 = p(t(:,3),1); y3 = p(t(:,3),2);
twiceArea = (x2-x1).*(y3-y1)-(y2-y1).*(x3-x1);
if any(twiceArea<=0)
    error('Triangles must be counterclockwise and nondegenerate.');
end
area = twiceArea/2;
b = [y2-y3,y3-y1,y1-y2];
c = [x3-x2,x1-x3,x2-x1];

I = reshape(repelem(t,1,3).',[],1);
J = reshape(repmat(t,1,3).',[],1);
Kv = zeros(9*nElem,1);
Mv = zeros(9*nElem,1);
massPattern = [2 1 1;1 2 1;1 1 2];
for a = 1:3
    for d = 1:3
        ids = (a-1)*3+d:9:9*nElem;
        Kv(ids) = (b(:,a).*b(:,d)+c(:,a).*c(:,d))./(4*area);
        Mv(ids) = area*(massPattern(a,d)/12);
    end
end

Kfull = sparse(I,J,Kv,nNode,nNode);
Mfull = sparse(I,J,Mv,nNode,nNode);
Dfull = accumarray(t(:),repmat(area/3,3,1),[nNode,1]);
K = Kfull(free,free);
M = Mfull(free,free);
D = Dfull(free);

fem = struct('K',(K+K')/2,'M',(M+M')/2,'D',D, ...
    'Mlump',D,'Kfull',(Kfull+Kfull')/2,'Mfull',(Mfull+Mfull')/2, ...
    'Dfull',Dfull,'MlumpFull',Dfull,'area',area);

fem.trueLumpIdentityError = norm(Dfull-Mfull*ones(nNode,1))/norm(Dfull);
fem.reducedRowSumRelativeDifference = norm(D-M*ones(numel(free),1))/norm(D);
end
