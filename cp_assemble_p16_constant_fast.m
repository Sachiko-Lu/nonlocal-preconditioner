function [A,stats] = cp_assemble_p16_constant_fast(mesh,delta,geomTol)
if nargin<3,geomTol=2e-11;end
validateattributes(delta,{'double'},{'scalar','positive','finite'});
validateattributes(geomTol,{'double'},{'scalar','positive','finite'});
tStart=tic;
p=mesh.p;t=mesh.t;ne=size(t,1);nd=numel(mesh.free);
assert(nd<=20000,'This implementation uses a dense coupling accumulator.');
px=reshape(p(t',1),3,ne)';py=reshape(p(t',2),3,ne)';
area=abs((px(:,2)-px(:,1)).*(py(:,3)-py(:,1))- ...
    (px(:,3)-px(:,1)).*(py(:,2)-py(:,1)))/2;
assert(all(area>0));
dofs=reshape(mesh.full2free(t),size(t));
theta=(2*(0:15)'+1)*pi/16;N=[cos(theta),sin(theta)];
support=delta*cos(pi/16);polygonArea=8*sin(pi/8);
lo=zeros(ne,16);hi=lo;
for k=1:16
    projected=px*N(k,1)+py*N(k,2);
    lo(:,k)=min(projected,[],2);hi(:,k)=max(projected,[],2);
end
bbox=[min(px,[],2),max(px,[],2),min(py,[],2),max(py,[],2)];
Im=zeros(9*ne,1);Jm=Im;Vm=Im;
for a=1:3
    for b=1:3
        ix=((a-1)*3+b-1)*ne+(1:ne);
        Im(ix)=dofs(:,a);Jm(ix)=dofs(:,b);
        Vm(ix)=area*(1+(a==b))/12;
    end
end
keep=Im>0 & Jm>0;M=sparse(Im(keep),Jm(keep),Vm(keep),nd,nd);
clear Im Jm Vm keep
C=zeros(nd,nd);
stats=struct('nCandidatePairs',0,'nFullPairs',0,'nCutPairs',0, ...
    'nNonemptyCutPairs',0,'nBlocks',0,'nFourSimplices',0, ...
    'nQuadPoints',0,'maxVolumeRelativeError',0, ...
    'maxConstraintViolation',0,'maxBarycentricViolation',0, ...
    'minWeight',inf,'interiorPairVolume',0,'geometryTolerance',geomTol, ...
    'nDelaunayFallbacks',0);
absoluteTol=64*eps(max(1,max(abs(p(:)))+delta));
for ie=1:ne
    idsT=dofs(ie,:);activeT=idsT>0;
    cand=find((1:ne)'>=ie & bbox(:,2)>=bbox(ie,1)-delta & ...
        bbox(:,1)<=bbox(ie,2)+delta & bbox(:,4)>=bbox(ie,3)-delta & ...
        bbox(:,3)<=bbox(ie,4)+delta);
    stats.nCandidatePairs=stats.nCandidatePairs+numel(cand);
    if isempty(cand),continue;end
    minimum=lo(cand,:)-hi(ie,:);
    maximum=hi(cand,:)-lo(ie,:);
    outside=any(minimum>support+absoluteTol,2);
    full=all(maximum<=support,2) & ~outside;
    jf=cand(full);
    stats.nFullPairs=stats.nFullPairs+numel(jf);
    if ~isempty(jf)
        multiplicity=2-((jf==ie));
        stats.interiorPairVolume=stats.interiorPairVolume+ ...
            area(ie)*sum(multiplicity.*area(jf));
        self=jf==ie;
        if any(self) && any(activeT)
            dt=idsT(activeT);
            C(dt,dt)=C(dt,dt)+(area(ie)^2/9)*ones(numel(dt),numel(dt));
            stats.nBlocks=stats.nBlocks+1;
        end
        other=jf(~self);
        if ~isempty(other) && any(activeT)
            ds=dofs(other,:);weights=repmat(area(ie)*area(other)/9,1,3);
            valid=ds>0;
            fullIds=ds(valid);fullWeights=weights(valid);
            s=accumarray(fullIds(:),fullWeights(:),[nd,1]);
            jj=find(s);values=s(jj);
            for a=find(activeT)
                ii=idsT(a);
                C(ii,jj)=C(ii,jj)+values';
                C(jj,ii)=C(jj,ii)+values;
            end
            stats.nBlocks=stats.nBlocks+sum(any(ds>0,2));
        end
    end
    cut=cand(~outside & ~full);
    stats.nCutPairs=stats.nCutPairs+numel(cut);
    T=[px(ie,:)',py(ie,:)'];
    for je=reshape(cut,1,[])
        S=[px(je,:)',py(je,:)'];
        [Q,volume,qa]=clipped_coupling(T,S,N,support,geomTol);
        if volume==0,continue;end
        stats.nNonemptyCutPairs=stats.nNonemptyCutPairs+1;
        stats.interiorPairVolume=stats.interiorPairVolume+ ...
            (2-(je==ie))*volume;
        stats.nFourSimplices=stats.nFourSimplices+qa.nSimplices;
        stats.nQuadPoints=stats.nQuadPoints+5*qa.nSimplices;
        stats.maxVolumeRelativeError=max(stats.maxVolumeRelativeError,qa.volumeError);
        stats.maxConstraintViolation=max(stats.maxConstraintViolation,qa.constraintError);
        stats.maxBarycentricViolation=max(stats.maxBarycentricViolation,qa.barycentricError);
        stats.minWeight=min(stats.minWeight,qa.minWeight);
        stats.nDelaunayFallbacks=stats.nDelaunayFallbacks+qa.usedDelaunayFallback;
        idsS=dofs(je,:);activeS=idsS>0;
        if ~any(activeT) || ~any(activeS),continue;end
        dt=idsT(activeT);ds=idsS(activeS);block=Q(activeT,activeS);
        if ie==je
            block=tril(block)+tril(block,-1)';
            C(dt,ds)=C(dt,ds)+block;
        else
            C(dt,ds)=C(dt,ds)+block;
            C(ds,dt)=C(ds,dt)+block';
        end
        stats.nBlocks=stats.nBlocks+1;
    end
    if mod(ie,max(1,round(ne/20)))==0 || ie==ne
        fprintf('  P16 coupling %d/%d elements: full=%d, cut=%d, %.1fs\n', ...
            ie,ne,stats.nFullPairs,stats.nCutPairs,toc(tStart));
    end
end
stats.couplingSymmetry=norm(C-C','fro')/max(norm(C,'fro'),realmin);
assert(stats.couplingSymmetry<1e-12,'Coupling matrix lost symmetry.');
C=tril(C)+tril(C,-1)';
A=2*(4/pi)*delta^-4*(delta^2*polygonArea*M-sparse(C));
clear C
stats.seconds=toc(tStart);stats.ndof=nd;stats.nnzA=nnz(A);
stats.symmetry=norm(A-A','fro')/max(norm(A,'fro'),realmin);
stats.domainArea=sum(area);stats.polygonArea=polygonArea;
stats.pairVolume=2*sum(area)*delta^2*polygonArea-stats.interiorPairVolume;
stats.exactPairVolume=NaN;stats.relativePairVolumeError=NaN;
if abs(sum(area)-1)<1e-10 && all(p(:)>=-1e-12 & p(:)<=1+1e-12) ...
        && exist('cp_p16_pair_volume','file')==2
    stats.exactPairVolume=cp_p16_pair_volume(delta);
    stats.relativePairVolumeError=abs(stats.pairVolume-stats.exactPairVolume) ...
        /stats.exactPairVolume;
end
if isinf(stats.minWeight),stats.minWeight=NaN;end
stats.assemblyMethod='exact mass minus interior coupling; positive degree-two clipped-pair cubature';
end

function [Q,volume,qa]=clipped_coupling(T,S,N,support,tol)
Q=zeros(3);volume=0;
qa=struct('nSimplices',0,'volumeError',0,'constraintError',0, ...
    'barycentricError',0,'minWeight',inf,'usedDelaunayFallback',false);
ct=mean(T,1);cs=mean(S,1);
scale=max([vecnorm(T([2 3 1],:)-T,2,2);vecnorm(S([2 3 1],:)-S,2,2)]);
T=(T-ct)/scale;S=(S-cs)/scale;
rhs=(support-N*(cs-ct)')/scale;
lower=min(S*N',[],1)'-max(T*N',[],1)';
upper=max(S*N',[],1)'-min(T*N',[],1)';
if any(lower>rhs+tol),return;end
active=upper>rhs-tol;
[AT,bT]=cp_triangle_halfspace(T);[AS,bS]=cp_triangle_halfspace(S);
H=[AT zeros(3,2);zeros(3,2) AS;-N(active,:) N(active,:)];
b=[bT;bS;rhs(active)];
norms=vecnorm(H,2,2);H=H./norms;b=b./norms;
comb=constraint_combinations(size(H,1));
V=zeros(size(comb,1),4);nv=0;
for j=1:size(comb,1)
    Hsub=H(comb(j,:),:);
    if rcond(Hsub)<1e-12,continue;end
    v=Hsub\b(comb(j,:));
    if all(H*v<=b+tol),nv=nv+1;V(nv,:)=v';end
end
if nv<5,return;end
V=uniquetol(V(1:nv,:),10*tol,'ByRows',true,'DataScale',1);
if size(V,1)<5 || rank(V-mean(V,1),2e-11)<4,return;end
try
    tr=cp_triangulate_scaled4d(V, ...
        struct('rankTolerance',2e-12,'volumeTolerance',2e-8));
catch failure
    failureDir=fullfile(fileparts(mfilename('fullpath')),'geometry_failures');
    if ~isfolder(failureDir),mkdir(failureDir);end
    failurePath=[tempname(failureDir),'.mat'];
    save(failurePath,'V','H','b','T','S','N','support','tol','scale','ct','cs','rhs');
    warning('Saved the failed geometry to %s.',failurePath);
    rethrow(failure)
end
if ~tr.isFullDimensional,return;end
aa=(1+4/sqrt(6))/5;bb=(1-1/sqrt(6))/5;
for j=1:numel(tr.simplexVolumes)
    R=squeeze(tr.simplexVertices(j,:,:));
    xy=(aa-bb)*R+bb*sum(R,1);
    lt=cp_barycentric(T,xy(:,1:2));ls=cp_barycentric(S,xy(:,3:4));
    weight=tr.simplexVolumes(j)/5;
    Q=Q+weight*(lt'*ls);
    qa.barycentricError=max(qa.barycentricError, ...
        max([0;-lt(:);-ls(:);lt(:)-1;ls(:)-1]));
end
Q=scale^4*Q;volume=scale^4*tr.volume;
qa.nSimplices=numel(tr.simplexVolumes);qa.volumeError=tr.volumeRelativeError;
qa.minWeight=scale^4*min(tr.simplexVolumes)/5;
qa.constraintError=max(0,max(max(H*V'-b)));
qa.usedDelaunayFallback=tr.usedDelaunayFallback;
end

function comb=constraint_combinations(n)
persistent saved
if isempty(saved),saved=cell(22,1);end
if isempty(saved{n})
    c=nchoosek(1:n,4);
    keep=sum(c<=3,2)<=2 & sum(c>3 & c<=6,2)<=2 & sum(c>6,2)<=2;
    saved{n}=c(keep,:);
end
comb=saved{n};
end
