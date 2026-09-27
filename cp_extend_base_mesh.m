function mesh = cp_extend_base_mesh(base,deltaMax,nExterior)
arguments
    base (1,1) struct
    deltaMax (1,1) double {mustBePositive}
    nExterior (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(nExterior,2)}
end

if ~isfield(base,'p') || ~isfield(base,'t') || ~isfield(base,'free')
    error('base must contain p, t, and free.');
end

aux=cp_make_extended_mesh(nExterior,deltaMax,314159);
extIds=find(~aux.isInteriorElem);
used=unique(aux.t(extIds,:));
old2new=zeros(size(aux.p,1),1);
old2new(used)=size(base.p,1)+(1:numel(used));

p=[base.p;aux.p(used,:)];
tInterior=base.t;
tExterior=reshape(old2new(aux.t(extIds,:)),size(aux.t(extIds,:)));
t=[tInterior;tExterior];
isInteriorElem=[true(size(tInterior,1),1);false(size(tExterior,1),1)];

free=base.free(:);
full2free=zeros(size(p,1),1);
full2free(free)=1:numel(free);

x1=p(t(:,1),1);x2=p(t(:,2),1);x3=p(t(:,3),1);
y1=p(t(:,1),2);y2=p(t(:,2),2);y3=p(t(:,3),2);
bbox=[min([x1,x2,x3],[],2),max([x1,x2,x3],[],2), ...
      min([y1,y2,y3],[],2),max([y1,y2,y3],[],2)];
edgeLengths=[hypot(x1-x2,y1-y2);hypot(x2-x3,y2-y3);hypot(x3-x1,y3-y1)];

if isfield(base,'boundary')
    boundary=base.boundary(:);
elseif isfield(base,'isBoundary')
    boundary=find(base.isBoundary);
else
    boundary=setdiff((1:size(base.p,1))',free);
end

mesh=struct('p',p,'t',t,'isInteriorElem',isInteriorElem,'free',free, ...
    'boundary',boundary,'full2free',full2free,'bbox',bbox, ...
    'h',1/nExterior,'hmax',max(edgeLengths),'deltaMax',deltaMax, ...
    'nInteriorElements',size(tInterior,1),'nExteriorElements',size(tExterior,1));
end
