function part = cp_kernel_partition(m,q,ker,opts)
arguments
    m (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(m,4)}
    q (1,1) double {mustBePositive}
    ker (1,1) struct
    opts.massRelTol (1,1) double {mustBePositive} = 2e-5
    opts.maxMassBoxes (1,1) double {mustBeInteger,mustBePositive} = 50000
end
if mod(m,2)~=0, error('m must be even.'); end
theta=((0:m-1)'+0.5)*(2*pi/m);
v=[cos(theta),sin(theta)];
Z=zeros(m,3,2);
for j=1:m
    jp=mod(j,m)+1;
    Z(j,:,:)=[0 0;v(j,:);v(jp,:)];
end
while max_triangle_diameter(Z)>q*(1+64*eps)
    Z=red_refine_triangles(Z);
end
nz=size(Z,1);
cent=squeeze(mean(Z,2));
mu=zeros(nz,1); lo=mu; hi=mu; nbox=zeros(nz,1); partner=zeros(nz,1);
done=false(nz,1);
for j=1:nz
    if done(j), continue; end
    [~,k]=min(vecnorm(cent+cent(j,:),2,2));
    if norm(cent(k,:)+cent(j,:))>1e-9
        error('Could not pair centrally symmetric kernel cells.');
    end
    tri=squeeze(Z(j,:,:));
    [mu(j),lo(j),hi(j),nbox(j)]=double_range_mass(tri,ker,opts.massRelTol,opts.maxMassBoxes);
    mu(k)=mu(j); lo(k)=lo(j); hi(k)=hi(j); nbox(k)=nbox(j);
    partner([j k])=[k j]; done([j k])=true;
end
area=triangle_area_batch(Z);
part=struct('m',m,'qTarget',q,'qActual',max_triangle_diameter(Z), ...
    'triangles',Z,'centres',cent,'area',area,'mu',mu,'muLo',lo,'muHi',hi, ...
    'rhoAvg',mu./area,'partner',partner,'massBoxes',nbox, ...
    'massSum',sum(mu),'massRelRadius',max((hi-lo)./(2*max(lo,realmin))), ...
    'kernel',ker.id,'beta',ker.beta,'kernelMass',ker.mass, ...
    'kernelRstar',ker.rStar);
end

function Z2=red_refine_triangles(Z)
n=size(Z,1); Z2=zeros(4*n,3,2);
for k=1:n
    T=squeeze(Z(k,:,:));
    a=T(1,:);b=T(2,:);c=T(3,:);
    ab=(a+b)/2;bc=(b+c)/2;ca=(c+a)/2;
    Z2(4*k-3,:,:)= [a;ab;ca];
    Z2(4*k-2,:,:)= [ab;b;bc];
    Z2(4*k-1,:,:)= [ca;bc;c];
    Z2(4*k,:,:)=   [ab;bc;ca];
end
end

function d=max_triangle_diameter(Z)
a=squeeze(Z(:,1,:));b=squeeze(Z(:,2,:));c=squeeze(Z(:,3,:));
d=max([vecnorm(a-b,2,2);vecnorm(b-c,2,2);vecnorm(c-a,2,2)]);
end

function a=triangle_area_batch(Z)
p=squeeze(Z(:,1,:));q=squeeze(Z(:,2,:));r=squeeze(Z(:,3,:));
a=abs((q(:,1)-p(:,1)).*(r(:,2)-p(:,2))-(q(:,2)-p(:,2)).*(r(:,1)-p(:,1)))/2;
end

function [mid,lo,hi,nbox]=double_range_mass(T,ker,relTol,maxBoxes)
J=abs(det([T(2,:)-T(1,:);T(3,:)-T(1,:)]));
if strcmp(ker.id,'K1')
    mid=(4/pi)*(J/2); lo=mid; hi=mid; nbox=1; return
end
boxes=[0 1 0 1];
[blo,bhi]=box_bounds(boxes(1,:),T,ker,J);
los=blo; his=bhi;
while true
    lo=sum(los); hi=sum(his); mid=(lo+hi)/2;
    if (hi-lo)<=2*relTol*max(mid,1e-15) || size(boxes,1)>=maxBoxes
        break
    end
    [~,k]=max(his-los);
    B=boxes(k,:); ds=B(2)-B(1); dt=B(4)-B(3);
    if ds>=dt
        z=(B(1)+B(2))/2;
        B1=[B(1) z B(3) B(4)]; B2=[z B(2) B(3) B(4)];
    else
        z=(B(3)+B(4))/2;
        B1=[B(1) B(2) B(3) z]; B2=[B(1) B(2) z B(4)];
    end
    [l1,h1]=box_bounds(B1,T,ker,J); [l2,h2]=box_bounds(B2,T,ker,J);
    boxes(k,:)=B1; los(k)=l1; his(k)=h1;
    boxes(end+1,:)=B2; los(end+1,1)=l2; his(end+1,1)=h2; %#ok<AGROW>
end
nbox=size(boxes,1);
inflate=256*eps*max(1,abs(mid));
lo=max(0,lo-inflate); hi=hi+inflate;
midGauss=duffy_gauss_mass(T,ker,14);
mid=min(hi,max(lo,midGauss));
end

function [L,H]=box_bounds(B,T,ker,J)
s0=B(1);s1=B(2);t0=B(3);t1=B(4);
a=T(1,:); u=T(2,:)-a; v=T(3,:)-a;
[S,Tp]=ndgrid([s0 s1],[t0 t1]);
pts=zeros(4,2); n=0;
for i=1:2
    for j=1:2
        n=n+1; pts(n,:)=a+S(i,j)*((1-Tp(i,j))*u+Tp(i,j)*v);
    end
end
xlo=min(pts(:,1));xhi=max(pts(:,1)); ylo=min(pts(:,2));yhi=max(pts(:,2));
dx=max([xlo,0,-xhi]); dy=max([ylo,0,-yhi]);
rlo=hypot(dx,dy);
rhi=max(hypot(pts(:,1),pts(:,2)));
area=(s1-s0)*(t1-t0);
switch ker.id
    case 'K2'
        c=20/pi;
        rhoLo=c*max(0,1-rhi); rhoHi=c*max(0,1-rlo);
        gLo=J*s0*rhoLo; gHi=J*s1*rhoHi;
    case 'K3'
        c=3/pi;
        if norm(a)<1e-14
            wpts=[(1-t0)*T(2,:)+t0*T(3,:);(1-t1)*T(2,:)+t1*T(3,:)];
            wxlo=min(wpts(:,1));wxhi=max(wpts(:,1));wylo=min(wpts(:,2));wyhi=max(wpts(:,2));
            ddx=max([wxlo,0,-wxhi]); ddy=max([wylo,0,-wyhi]);
            wrlo=max(hypot(ddx,ddy),realmin); wrhi=max(hypot(wpts(:,1),wpts(:,2)));
            gLo=c*J/wrhi; gHi=c*J/wrlo;
        else
            rlo=max(rlo,realmin);
            gLo=c*J*s0/rhi; gHi=c*J*s1/rlo;
        end
    otherwise
        error('Validated enclosure implemented only for K1--K3.');
end
L=max(0,area*gLo); H=max(L,area*gHi);
pad=128*eps*max(1,H); L=max(0,L-pad); H=H+pad;
end

function val=duffy_gauss_mass(T,ker,n)
[x,w]=gauss01(n);
[S,Q]=ndgrid(x,x); [WS,WQ]=ndgrid(w,w);
a=T(1,:);u=T(2,:)-a;v=T(3,:)-a;
xx=a(1)+S.*((1-Q)*u(1)+Q*v(1));
yy=a(2)+S.*((1-Q)*u(2)+Q*v(2));
J=abs(det([u;v]));
r=hypot(xx,yy);
switch ker.id
    case 'K1', rho=(4/pi)*ones(size(r));
    case 'K2', rho=(20/pi)*max(0,1-r);
    case 'K3', rho=(3/pi)./r;
end
val=sum(WS.*WQ.*rho.*S*J,'all');
end

function [x,w]=gauss01(n)
k=(1:n-1)'; b=k./sqrt(4*k.^2-1);
[V,D]=eig(diag(b,1)+diag(b,-1));
[x,ix]=sort(diag(D)); V=V(:,ix); w=2*(V(1,:)'.^2);
x=(x+1)/2; w=w/2;
end
