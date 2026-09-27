function part=cp_p16_partition(q,alpha)
if nargin<2, alpha=0; end
assert(alpha>=0 && alpha<1);
part=cp_kernel_partition(16,q,cp_kernel('K1'));
a=-pi/16; R=[cos(a),-sin(a);sin(a),cos(a)];
for j=1:size(part.triangles,1)
    T=squeeze(part.triangles(j,:,:))*R';
    T(abs(T)<64*eps)=0;
    part.triangles(j,:,:)=T;
end
part.centres=squeeze(mean(part.triangles,2));
sector=sign(part.centres(:,1).^2-part.centres(:,2).^2);
assert(all(abs(sector)==1));
moment=zeros(2);
for j=1:numel(sector)
    T=squeeze(part.triangles(j,:,:));
    signs=T(:,1).^2-T(:,2).^2;
    assert(all(sector(j)*signs>=-1e-13), ...
        'A kernel triangle crosses an angular-kernel jump line.');
    part.area(j)=abs(det([T(2,:)-T(1,:);T(3,:)-T(1,:)]))/2;
    density=(4/pi)*(1+alpha*sector(j));
    s=sum(T,1);
    moment=moment+density*part.area(j)/12*(s'*s+T'*T);
end
part.rhoAvg=(4/pi)*(1+alpha*sector);
part.mu=part.area.*part.rhoAvg;
part.muLo=part.mu; part.muHi=part.mu;
part.massSum=sum(part.mu); part.massRelRadius=0;
part.massBoxes=ones(numel(part.mu),1);
part.kernel=sprintf('angular_a%.1f',alpha);
if alpha==0, part.kernel='K1'; end
part.beta=1/8; part.kernelMass=4; part.kernelRstar=1;
theta=(0:15)'*pi/8;
part.vertices=[cos(theta),sin(theta)];
part.vertices(abs(part.vertices)<64*eps)=0;
part.polygonArea=8*sin(pi/8);
part.secondMomentTensor=moment;
part.localLimitCoefficient=moment(1,1);
part.alpha=alpha;
assert(abs(part.massSum-(4/pi)*part.polygonArea)<1e-12);
assert(max(abs(part.rhoAvg-part.rhoAvg(part.partner)))<1e-14);
end
