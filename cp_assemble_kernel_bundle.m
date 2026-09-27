function bundle = cp_assemble_kernel_bundle(mesh,partitions,delta)
t0=tic;
nKernel=numel(partitions); part=partitions{1};
for k=2:nKernel
    assert(isequal(part.triangles,partitions{k}.triangles));
end
assert(delta<=mesh.deltaMax*(1+1e-12));
p=mesh.p; t=mesh.t; ne=size(t,1); nd=numel(mesh.free);
nz=size(part.triangles,1); zb=zeros(nz,4);
rho=zeros(nz,nKernel);
for iz=1:nz
    Z=squeeze(part.triangles(iz,:,:));
    zb(iz,:)=[min(Z(:,1)),max(Z(:,1)),min(Z(:,2)),max(Z(:,2))];
end
for k=1:nKernel, rho(:,k)=partitions{k}.rhoAvg; end
assert(all(rho(:)>0));
Icell=cell(10000,1); Jcell=Icell; Vcell=cell(10000,nKernel); nBlocks=0;
stats=struct('candidateElementPairs',0,'candidateTriples',0, ...
    'nonemptyPairCells',0,'nFourSimplices',0,'nQuadPoints',0, ...
    'pairVolume',0,'minWeight',inf,'maxConstraintViolation',0, ...
    'maxBarycentricViolation',0,'skippedZeroDof',0);
for ie=1:ne
    bt=mesh.bbox(ie,:); sb=mesh.bbox;
    cand=(sb(:,2)-bt(1)>=-delta)&(sb(:,1)-bt(2)<=delta)& ...
        (sb(:,4)-bt(3)>=-delta)&(sb(:,3)-bt(4)<=delta);
    if ~mesh.isInteriorElem(ie), cand=cand&mesh.isInteriorElem; end
    T=p(t(ie,:),:);
    for je=reshape(find(cand),1,[])
        if ~(mesh.isInteriorElem(ie)||mesh.isInteriorElem(je)), continue; end
        stats.candidateElementPairs=stats.candidateElementPairs+1;
        bs=mesh.bbox(je,:);
        db=[(bs(1)-bt(2))/delta,(bs(2)-bt(1))/delta, ...
            (bs(3)-bt(4))/delta,(bs(4)-bt(3))/delta];
        zcand=find(zb(:,2)>=db(1)&zb(:,1)<=db(2)& ...
            zb(:,4)>=db(3)&zb(:,3)<=db(4));
        S=p(t(je,:),:);
        for iz=reshape(zcand,1,[])
            stats.candidateTriples=stats.candidateTriples+1;
            Z=squeeze(part.triangles(iz,:,:));
            pc=cp_paircell_4d(T,S,Z,delta,tol=5e-11);
            if pc.empty, continue; end
            stats.nonemptyPairCells=stats.nonemptyPairCells+1;
            stats.nFourSimplices=stats.nFourSimplices+pc.nSimplices;
            stats.nQuadPoints=stats.nQuadPoints+numel(pc.w);
            stats.pairVolume=stats.pairVolume+pc.volume;
            stats.minWeight=min(stats.minWeight,min(pc.w));
            stats.maxConstraintViolation=max(stats.maxConstraintViolation, ...
                max(max(pc.A*pc.q'-pc.b)));
            lamT=cp_barycentric(T,pc.q(:,1:2));
            lamS=cp_barycentric(S,pc.q(:,3:4));
            stats.maxBarycentricViolation=max(stats.maxBarycentricViolation, ...
                max([0;-lamT(:);-lamS(:);lamT(:)-1;lamS(:)-1]));
            idsT=mesh.full2free(t(ie,:)); idsS=mesh.full2free(t(je,:));
            dofs=unique([reshape(idsT(idsT>0),[],1);reshape(idsS(idsS>0),[],1)]);
            dofs=dofs(:)';
            if isempty(dofs)
                stats.skippedZeroDof=stats.skippedZeroDof+1; continue
            end
            D=zeros(numel(pc.w),numel(dofs));
            for a=1:3
                if idsT(a)>0, D(:,dofs==idsT(a))=D(:,dofs==idsT(a))-lamT(:,a); end
                if idsS(a)>0, D(:,dofs==idsS(a))=D(:,dofs==idsS(a))+lamS(:,a); end
            end
            R=delta^-4*(D'*(pc.w.*D));
            R=tril(R)+tril(R,-1)';
            if norm(R,'fro')<=1e-30, continue; end
            nBlocks=nBlocks+1;
            if nBlocks>numel(Icell)
                Icell=[Icell;cell(10000,1)]; %#ok<AGROW>
                Jcell=[Jcell;cell(10000,1)]; %#ok<AGROW>
                Vcell=[Vcell;cell(10000,nKernel)]; %#ok<AGROW>
            end
            [ii,jj]=ndgrid(dofs,dofs);
            Icell{nBlocks}=ii(:); Jcell{nBlocks}=jj(:);
            for k=1:nKernel, Vcell{nBlocks,k}=rho(iz,k)*R(:); end
        end
    end
    if mod(ie,max(1,round(ne/10)))==0
        fprintf('  shared geometry %d/%d elements, %d blocks\n',ie,ne,nBlocks);
    end
end
ii=vertcat(Icell{1:nBlocks}); jj=vertcat(Jcell{1:nBlocks});
bundle=cell(nKernel,1);
for k=1:nKernel
    A=sparse(ii,jj,vertcat(Vcell{1:nBlocks,k}),nd,nd);
    st=stats; st.nBlocks=nBlocks; st.ndof=nd; st.nnzA=nnz(A);
    st.symmetry=norm(A-A','fro')/max(norm(A,'fro'),realmin);
    bundle{k}=struct('Adbl',A,'stats',st,'delta',delta, ...
        'beta',partitions{k}.beta,'part',partitions{k}, ...
        'kernel',partitions{k}.kernel,'m',part.m,'q',part.qActual, ...
        'geometryCertified',false,'massEnclosureIsFormal',false);
end
for k=1:nKernel, bundle{k}.stats.seconds=toc(t0); end
end
