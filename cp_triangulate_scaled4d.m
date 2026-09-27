function out=cp_triangulate_scaled4d(vertices,opts)
if nargin<2,opts=struct;end
center=mean(vertices,1);
[~,S,Q]=svd(vertices-center,0);
scales=diag(S);
assert(numel(scales)==4 && all(scales>0),'The vertex set must have affine rank four.');
normalized=((vertices-center)*Q)./scales';
try
    out=geometry_triangulate_polytope4d(normalized,opts);
    out.usedDelaunayFallback=false;
catch
    try
        simplices=delaunayn(normalized,{'Qbb','Qc','Qz','Qt'});
    catch
        simplices=delaunayn(normalized,{'QJ','Qbb','Qc'});
    end
    try
        [facets,hullVolume]=convhulln(normalized,{'Qt','Qx'});
    catch
        [facets,hullVolume]=convhulln(normalized,{'QJ'});
    end
    ns=size(simplices,1);simplexVertices=zeros(ns,5,4);volumes=zeros(ns,1);
    for j=1:ns
        R=normalized(simplices(j,:),:);
        simplexVertices(j,:,:)=R;
        volumes(j)=abs(det(R(2:5,:)-R(1,:)))/24;
    end
    keep=volumes>0;
    simplexVertices=simplexVertices(keep,:,:);volumes=volumes(keep);
    volumeError=abs(sum(volumes)-hullVolume)/max(hullVolume,realmin);
    if isfield(opts,'volumeTolerance'),volumeTolerance=opts.volumeTolerance;
    else,volumeTolerance=2e-8;end
    assert(volumeError<=volumeTolerance,'Scaled Delaunay partition failed the volume check.');
    out=struct('vertices',normalized,'facets',facets,'center',mean(normalized,1), ...
        'simplexVertices',simplexVertices,'simplexVolumes',volumes, ...
        'volume',sum(volumes),'hullVolume',hullVolume,'volumeRelativeError',volumeError, ...
        'isFullDimensional',true,'usedDelaunayFallback',true);
end
if ~out.isFullDimensional,return;end
jacobian=prod(scales);
for j=1:numel(out.simplexVolumes)
    R=squeeze(out.simplexVertices(j,:,:));
    out.simplexVertices(j,:,:)=(R.*scales')*Q'+center;
end
out.simplexVolumes=jacobian*out.simplexVolumes;
out.volume=jacobian*out.volume;
out.hullVolume=jacobian*out.hullVolume;
out.vertices=vertices;
out.center=center;
out.affineConditionNumber=max(scales)/min(scales);
end
