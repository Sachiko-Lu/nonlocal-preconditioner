function results=run_convergence_experiment(nValues,outputDir,saveCheckpoints)
if nargin<1,nValues=[4 8 16 32 64];end
if nargin<2,outputDir=fullfile(fileparts(mfilename('fullpath')),'results','section_6_1');end
if nargin<3,saveCheckpoints=false;end
if ~isfolder(outputDir),mkdir(outputDir);end
mFactor=4;geomTol=2e-11;
delta=0.2; rows=struct([]);
sourceOrder=20;sourceCheckOrder=32;
sourceOuterOrder=24;sourceOuterCheckOrder=32;
sourcePointCheck=verify_forcing(delta,sourceOrder,sourceCheckOrder);
fprintf('Homogeneous disk source check: %.3e\n',sourcePointCheck);
for n=nValues
    m=mFactor*n;
    fprintf('\nDISK ACCURACY n=%d m=%d delta=%.8g\n',n,m,delta);
    mesh=uniform_extended_mesh(n,delta);
    t0=tic;
    [A,stats]=cached_assembly(mesh,n,m,delta,geomTol);
    [bsource,normu2]=source_load(mesh,delta,sourceOuterOrder,sourceOrder);
    assert(abs(normu2-1/900)<3e-15);
    b=bsource;
    [L,flag]=chol(A,'lower');assert(flag==0,'Matrix is not SPD.');
    U=L'\(L\b);
    [eL2,eInterp]=l2errors(mesh,U,normu2,6);
    [eL2check,~]=l2errors(mesh,U,normu2,10);
    [bsourceCheck,~]=source_load(mesh,delta,sourceOuterCheckOrder,sourceCheckOrder);
    UsourceCheck=L'\(L\bsourceCheck);
    [eL2SourceCheck,~]=l2errors(mesh,UsourceCheck,normu2,6);
    row=struct('n_subdivisions',n,'n_dof',numel(mesh.free), ...
        'hmax',sqrt(2)/n,'delta',delta,'m_polygon',m, ...
        'ball_area_defect',pi-0.5*m*sin(2*pi/m), ...
        'kernel_L1_error',(4/pi)*(pi-0.5*m*sin(2*pi/m)), ...
        'relative_L2_error',eL2,'observed_L2_rate',NaN, ...
        'relative_interpolation_L2_error',eInterp, ...
        'relative_residual',norm(A*U-b)/norm(b), ...
        'matrix_symmetry',norm(A-A','fro')/norm(A,'fro'), ...
        'minimum_cholesky_diagonal',min(diag(L)), ...
        'source_pointwise_check_difference',sourcePointCheck, ...
        'error_quadrature_difference',abs(eL2-eL2check), ...
        'source_quadrature_relative_difference',norm(bsource-bsourceCheck)/norm(bsource), ...
        'source_quadrature_error_difference',abs(eL2-eL2SourceCheck), ...
        'source_inner_quadrature_order',sourceOrder, ...
        'source_inner_check_order',sourceCheckOrder, ...
        'source_outer_quadrature_order',sourceOuterOrder, ...
        'source_outer_check_order',sourceOuterCheckOrder, ...
        'relative_pair_volume_error',stats.relativePairVolumeError, ...
        'maximum_simplex_volume_discrepancy',stats.maxHullVolumeError, ...
        'minimum_positive_weight_unit',stats.minimumWeight, ...
        'geometry_tolerance_unit',geomTol, ...
        'unique_element_pair_geometries',stats.nGeometry, ...
        'nonempty_unique_pairs',stats.nNonempty, ...
        'assembled_ordered_pair_count',stats.orderedPairs, ...
        'assembly_seconds',stats.seconds, ...
        'total_seconds',toc(t0));
    if ~isempty(rows)
        row.observed_L2_rate=log(rows(end).relative_L2_error/eL2)/log(rows(end).hmax/row.hmax);
    end
    rows=[rows;row]; %#ok<AGROW>
    results=struct2table(rows);
    writetable(results,fullfile(outputDir,'table_6_1.csv'));
    if saveCheckpoints
        save(fullfile(outputDir,sprintf('checkpoint_n%d_m%d.mat',n,m)), ...
            'A','b','U','row','stats','delta','n','m','mesh','-v7.3');
    end
    fprintf('RESULT N=%d L2=%.12e rate=%.6f residual=%.3e volumeQA=%.3e time=%.2fs\n', ...
        row.n_dof,eL2,row.observed_L2_rate,row.relative_residual, ...
        row.relative_pair_volume_error,row.total_seconds);
    assert(row.relative_residual<1e-10);
    assert(row.matrix_symmetry<1e-12);
    assert(row.relative_pair_volume_error<2e-7);
    assert(row.error_quadrature_difference<2e-12);
    assert(row.source_quadrature_relative_difference<1e-6);
    assert(row.source_quadrature_error_difference<1e-7);
end
disp(results(:,{'n_dof','m_polygon','hmax','relative_L2_error','observed_L2_rate'}));
end

function [A,stats]=cached_assembly(mesh,n,m,delta,tol)
t0=tic;nd=numel(mesh.free);h=1/n;r=delta/h;
Tbase={ [0 0;1 0;1 1], [0 0;1 1;0 1] };
angles=2*pi*(0:m-1)'/m;normals=[cos(angles),sin(angles)];
support=r*cos(pi/m);R=ceil(r)+1;
A=sparse(nd,nd);
[ix,iy]=ndgrid(0:n-1,0:n-1);ix=ix(:);iy=iy(:);
Ibuf={};Jbuf={};Vbuf={};bufferCount=0;
volume=0;orderedPairs=0;nGeometry=0;nNonempty=0;maxHullError=0;minWeight=inf;
for ot=1:2
    T=Tbase{ot};idsT=triangle_ids(ix,iy,ot,mesh.nExt,mesh.nx);
    for os=1:2
        S=Tbase{os};
        for dx=-R:R
            for dy=-R:R
                sx=ix+dx;sy=iy+dy;
                keep=sx>=-mesh.nExt & sx<n+mesh.nExt & sy>=-mesh.nExt & sy<n+mesh.nExt;
                if ~any(keep),continue;end
                nGeometry=nGeometry+1;
                [B,pcvol,qa]=pair_block(T,S,[dx dy],normals,support,tol);
                if pcvol==0,continue;end
                nNonempty=nNonempty+1;
                maxHullError=max(maxHullError,qa.volumeError);
                minWeight=min(minWeight,qa.minWeight);
                idsS=triangle_ids(sx(keep),sy(keep),os,mesh.nExt,mesh.nx);
                ids=[idsT(keep,:),idsS];dofs=mesh.full2free(ids);
                mult=1+double(sx(keep)<0 | sx(keep)>=n | sy(keep)<0 | sy(keep)>=n);
                orderedPairs=orderedPairs+sum(mult);
                volume=volume+pcvol*h^4*sum(mult);
                B=(4/pi)*r^-4*B;
                ic={};jc={};vc={};nb=0;
                for a=1:6
                    for bb=1:6
                        sel=dofs(:,a)>0 & dofs(:,bb)>0;
                        if ~any(sel) || B(a,bb)==0,continue;end
                        nb=nb+1;ic{nb}=dofs(sel,a);jc{nb}=dofs(sel,bb); %#ok<AGROW>
                        vc{nb}=B(a,bb)*mult(sel); %#ok<AGROW>
                    end
                end
                if nb>0
                    bufferCount=bufferCount+1;
                    Ibuf{bufferCount}=vertcat(ic{:}); %#ok<AGROW>
                    Jbuf{bufferCount}=vertcat(jc{:}); %#ok<AGROW>
                    Vbuf{bufferCount}=vertcat(vc{:}); %#ok<AGROW>
                end
                if bufferCount>=64
                    A=A+sparse(vertcat(Ibuf{:}),vertcat(Jbuf{:}),vertcat(Vbuf{:}),nd,nd);
                    Ibuf={};Jbuf={};Vbuf={};bufferCount=0;
                end
            end
            fprintf('  geometry ot=%d os=%d dx=%d/%d nonempty=%d time=%.1fs\n', ...
                ot,os,dx,R,nNonempty,toc(t0));
        end
    end
end
if bufferCount>0,A=A+sparse(vertcat(Ibuf{:}),vertcat(Jbuf{:}),vertcat(Vbuf{:}),nd,nd);end
exactVolume=ordered_pair_volume(m,delta);
stats=struct('nGeometry',nGeometry,'nNonempty',nNonempty,'orderedPairs',orderedPairs, ...
    'pairVolume',volume,'exactPairVolume',exactVolume, ...
    'relativePairVolumeError',abs(volume-exactVolume)/exactVolume, ...
    'maxHullVolumeError',maxHullError,'minimumWeight',minWeight,'seconds',toc(t0));
end

function [B,volume,qa]=pair_block(T,S,offset,N,support,tol)
B=zeros(6);volume=0;qa=struct('volumeError',0,'minWeight',inf);
lo=min(S*N',[],1)'-max(T*N',[],1)';
hi=max(S*N',[],1)'-min(T*N',[],1)';
rhs=support-N*offset';
if any(lo>=rhs+tol),return;end
active=hi>rhs-tol;
[AT,bT]=cp_triangle_halfspace(T);[AS,bS]=cp_triangle_halfspace(S);
A=[AT zeros(3,2);zeros(3,2) AS;-N(active,:) N(active,:)];
b=[bT;bS;rhs(active)];
nr=vecnorm(A,2,2);A=A./nr;b=b./nr;
comb=nchoosek(1:size(A,1),4);
valid=sum(comb<=3,2)<=2 & sum(comb>3 & comb<=6,2)<=2 & sum(comb>6,2)<=2;
comb=comb(valid,:);V=zeros(size(comb,1),4);nv=0;
for j=1:size(comb,1)
    C=A(comb(j,:),:);
    if rcond(C)<1e-12,continue;end
    v=C\b(comb(j,:));
    if all(A*v<=b+tol),nv=nv+1;V(nv,:)=v';end
end
V=V(1:nv,:);
if nv<5,return;end
V=uniquetol(V,10*tol,'ByRows',true,'DataScale',1);
if size(V,1)<5 || rank(V-mean(V,1),2e-11)<4,return;end
tr=geometry_triangulate_polytope4d(V,struct('rankTolerance',2e-12,'volumeTolerance',2e-8));
if ~tr.isFullDimensional,return;end
aa=(1+4/sqrt(6))/5;bb=(1-1/sqrt(6))/5;
for j=1:numel(tr.simplexVolumes)
    R=squeeze(tr.simplexVertices(j,:,:));
    Q=(aa-bb)*R+bb*sum(R,1);w=tr.simplexVolumes(j)/5;
    D=[-cp_barycentric(T,Q(:,1:2)),cp_barycentric(S,Q(:,3:4))];
    B=B+w*(D'*D);
end
B=tril(B)+tril(B,-1)';
volume=tr.volume;qa.volumeError=tr.volumeRelativeError;
qa.minWeight=min(tr.simplexVolumes)/5;
end

function ids=triangle_ids(i,j,o,nExt,nx)
a=(i+nExt+1)+nx.*(j+nExt);
if o==1,ids=[a a+1 a+nx+1];else,ids=[a a+nx+1 a+nx];end
end

function mesh=uniform_extended_mesh(n,delta)
h=1/n;nExt=ceil(delta/h)+1;xv=(-nExt:n+nExt)*h;
[X,Y]=ndgrid(xv,xv);p=[X(:),Y(:)];nx=numel(xv);
[i,j]=ndgrid(-nExt:n+nExt-1,-nExt:n+nExt-1);
t1=triangle_ids(i(:),j(:),1,nExt,nx);t2=triangle_ids(i(:),j(:),2,nExt,nx);
t=zeros(2*numel(i),3);t(1:2:end,:)=t1;t(2:2:end,:)=t2;
insideCell=i(:)>=0 & i(:)<n & j(:)>=0 & j(:)<n;
inside=repelem(insideCell,2);
free=find(all(p>100*eps & p<1-100*eps,2));
map=zeros(size(p,1),1);map(free)=1:numel(free);
px=reshape(p(t',1),3,[]);py=reshape(p(t',2),3,[]);
bbox=[min(px,[],1)' max(px,[],1)' min(py,[],1)' max(py,[],1)'];
mesh=struct('p',p,'t',t,'free',free,'full2free',map, ...
    'isInteriorElem',inside,'bbox',bbox,'hmax',sqrt(2)*h,'deltaMax',delta, ...
    'nExt',nExt,'nx',nx);
end

function [b,normu2]=source_load(mesh,delta,nq,sourceNq)
[lam,wref]=duffy_rule(nq);[sourceNodes,sourceWeights]=gauss01(sourceNq);
b=zeros(numel(mesh.free),1);normu2=0;
for e=reshape(find(mesh.isInteriorElem),1,[])
    ids=mesh.t(e,:);T=mesh.p(ids,:);xy=lam*T;
    w=abs(det([T(2,:)-T(1,:);T(3,:)-T(1,:)]))*wref;
    f=exactf(xy(:,1),xy(:,2),delta,sourceNodes,sourceWeights);
    normu2=normu2+sum(w.*exactu(xy(:,1),xy(:,2)).^2);
    dof=mesh.full2free(ids);
    for a=1:3
        if dof(a)>0,b(dof(a))=b(dof(a))+sum(w.*f.*lam(:,a));end
    end
end
end

function [eL2,eInterp]=l2errors(mesh,U,normu2,nq)
[lam,wref]=duffy_rule(nq);err=0;erri=0;
for e=reshape(find(mesh.isInteriorElem),1,[])
    ids=mesh.t(e,:);T=mesh.p(ids,:);xy=lam*T;
    w=abs(det([T(2,:)-T(1,:);T(3,:)-T(1,:)]))*wref;
    u=exactu(xy(:,1),xy(:,2));dof=mesh.full2free(ids);
    v=zeros(3,1);sel=dof>0;v(sel)=U(dof(sel));
    uh=lam*v;ui=lam*exactu(T(:,1),T(:,2));
    err=err+sum(w.*(u-uh).^2);erri=erri+sum(w.*(u-ui).^2);
end
eL2=sqrt(err/normu2);eInterp=sqrt(erri/normu2);
end

function volume=ordered_pair_volume(m,delta)
th=(2*(0:m-1)'+1)*pi/m;P=delta*[cos(th),sin(th)];
P=clip_positive(P,1);P=clip_positive(P,2);
[lam,wref]=duffy_rule(4);integral=0;c=mean(P,1);
for j=1:size(P,1)
    T=[c;P(j,:);P(mod(j,size(P,1))+1,:)];xy=lam*T;
    w=abs(det([T(2,:)-T(1,:);T(3,:)-T(1,:)]))*wref;
    integral=integral+sum(w.*(1-xy(:,1)).*(1-xy(:,2)));
end
volume=2*delta^2*(m/2*sin(2*pi/m))-4*integral;
end

function Q=clip_positive(P,dim)
Q=zeros(0,2);
for j=1:size(P,1)
    a=P(j,:);b=P(mod(j,size(P,1))+1,:);ina=a(dim)>=0;inb=b(dim)>=0;
    if ina,Q(end+1,:)=a;end %#ok<AGROW>
    if ina~=inb,Q(end+1,:)=a+(b-a)*(-a(dim)/(b(dim)-a(dim)));end %#ok<AGROW>
end
end

function discrepancy=verify_forcing(delta,nq,nqCheck)
[nodes,weights]=gauss01(nq);[nodesCheck,weightsCheck]=gauss01(nqCheck);
samples=[0.01 0.02;0.17 0.28;0.21 0.21;0.41 0.63;0.79 0.91;0.99 0.99];
f=exactf(samples(:,1),samples(:,2),delta,nodes,weights);
fCheck=exactf(samples(:,1),samples(:,2),delta,nodesCheck,weightsCheck);
fAdaptive=exactf_adaptive(samples(:,1),samples(:,2),delta);
interiorFormula=2*(profile(samples(4,1))+profile(samples(4,2)))-delta^2/3;
discrepancy=max([abs(f-fCheck);abs(f-fAdaptive);abs(f(4)-interiorFormula)]);
assert(discrepancy<1e-10);
end

function u=exactu(x,y)
u=x.*(1-x).*y.*(1-y);
end

function f=exactf(x,y,delta,nodes,weights)
shape=size(x);x=x(:);y=y(:);integral=zeros(size(x));
for k=1:numel(x)
    thetaLeft=asin(max(-1,-x(k)/delta));
    thetaRight=asin(min(1,(1-x(k))/delta));
    cuts=[thetaLeft;thetaRight];
    distances=[y(k),1-y(k)];
    for j=1:2
        if distances(j)<delta
            angle=acos(distances(j)/delta);
            candidates=[-angle;angle];
            cuts=[cuts;candidates(candidates>thetaLeft & candidates<thetaRight)]; %#ok<AGROW>
        end
    end
    cuts=unique(sort(cuts));
    for j=1:numel(cuts)-1
        a=cuts(j);c=cuts(j+1);
        theta=(c-a)*nodes+a;
        radial=delta*cos(theta);
        y1=x(k)+delta*sin(theta);
        lower=max(0,y(k)-radial);
        upper=min(1,y(k)+radial);
        values=profile(y1).*(profile_primitive(upper)-profile_primitive(lower)) ...
            .*delta.*cos(theta);
        integral(k)=integral(k)+(c-a)*sum(weights.*values);
    end
end
u=exactu(x,y);
f=8/(pi*delta^4)*(pi*delta^2*u-integral);
f=reshape(f,shape);
end

function f=exactf_adaptive(x,y,delta)
shape=size(x);x=x(:);y=y(:);integralValue=zeros(size(x));
for k=1:numel(x)
    left=max(0,x(k)-delta);right=min(1,x(k)+delta);
    cuts=[left;right];
    distances=[y(k),1-y(k)];
    for j=1:2
        if distances(j)<delta
            offset=sqrt(delta^2-distances(j)^2);
            candidates=[x(k)-offset;x(k)+offset];
            cuts=[cuts;candidates(candidates>left & candidates<right)]; %#ok<AGROW>
        end
    end
    cuts=unique(sort(cuts));
    for j=1:numel(cuts)-1
        integrand=@(y1) vertical_integrand(y1,x(k),y(k),delta);
        integralValue(k)=integralValue(k)+integral(integrand,cuts(j),cuts(j+1), ...
            'AbsTol',1e-14,'RelTol',1e-12);
    end
end
u=exactu(x,y);
f=8/(pi*delta^4)*(pi*delta^2*u-integralValue);
f=reshape(f,shape);
end

function value=vertical_integrand(y1,x,y,delta)
radial=sqrt(max(0,delta^2-(y1-x).^2));
lower=max(0,y-radial);upper=min(1,y+radial);
value=profile(y1).*(profile_primitive(upper)-profile_primitive(lower));
end

function value=profile(s)
value=s.*(1-s);
end

function value=profile_primitive(s)
value=s.^2/2-s.^3/3;
end

function [lambda,w]=duffy_rule(n)
[r,wr]=gauss01(n);[R,S]=ndgrid(r,r);[WR,WS]=ndgrid(wr,wr);
R=R(:);S=S(:);lambda=[1-R,R.*(1-S),R.*S];w=WR(:).*WS(:).*R;
end

function [x,w]=gauss01(n)
k=(1:n-1)';b=k./sqrt(4*k.^2-1);[V,D]=eig(diag(b,1)+diag(b,-1));
[x,ix]=sort(diag(D));V=V(:,ix);w=V(1,:)'.^2;x=(x+1)/2;
end
