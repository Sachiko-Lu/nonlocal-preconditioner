function [parts,kernelQA]=cp_p16_representative_kernels()

q=0.5;
template=cp_p16_partition(q,0);
stale={'secondMomentTensor','localLimitCoefficient','alpha'};
stale=intersect(stale,fieldnames(template));
if ~isempty(stale), template=rmfield(template,stale); end

names={"rho0","rhoS","rhoG","rho4","rhoC"};
descriptions={"constant radial", ...
    "discontinuous two-direction step", ...
    "smooth rotated mixed harmonic", ...
    "smooth rotated four-lobe", ...
    "aligned narrow double cone"};
nz=size(template.triangles,1);
area=template.area(:);
partner=template.partner(:);
c0=4/pi;
rhoStar=0.1*c0;
rho=zeros(nz,5);
rho(:,1)=c0;

sector=sign(template.centres(:,1).^2-template.centres(:,2).^2);
assert(all(abs(sector)==1),'A rhoS cell centre lies on a jump ray.');
for iz=1:nz
    T=squeeze(template.triangles(iz,:,:));
    jump=T(:,1).^2-T(:,2).^2;
    assert(all(sector(iz)*jump>=-2e-14), ...
        'A kernel cell crosses a rhoS jump ray.');
end
rho(:,2)=c0*(1+0.9*sector);

[bary,w]=positive_degree4_rule();
assert(all(w>0) && abs(sum(w)-1)<5e-15);
ruleError=degree4_rule_error(bary,w);
assert(ruleError<2e-14,'The triangle rule failed its degree-four audit.');
for iz=1:nz
    T=squeeze(template.triangles(iz,:,:));
    xy=bary*T;
    rho(iz,3)=sum(w.*kernel_values(3,xy));
    rho(iz,4)=sum(w.*kernel_values(4,xy));
end

R=19;
omega=pi/16;
phi=pi/16;
e=[cos(phi),sin(phi)];
eperp=[-sin(phi),cos(phi)];
coneScale=4/(pi+2*omega*(R-1));
cone=false(nz,1);
coneCrossings=0;
for iz=1:nz
    T=squeeze(template.triangles(iz,:,:));
    gv=tan(omega)*abs(T*e')-abs(T*eperp');
    gc=tan(omega)*abs(template.centres(iz,:)*e')- ...
        abs(template.centres(iz,:)*eperp');
    assert(abs(gc)>2e-14,'A rhoC cell centre lies on a cone boundary.');
    crosses=min(gv)<-2e-14 && max(gv)>2e-14;
    coneCrossings=coneCrossings+crosses;
    cone(iz)=gc>0;
    if cone(iz)
        assert(all(gv>=-2e-14),'A rhoC cell crosses a cone boundary.');
    else
        assert(all(gv<=2e-14),'A rhoC cell crosses a cone boundary.');
    end
end
assert(coneCrossings==0,'The rhoC cone is not aligned with the partition.');
assert(abs(sum(area(cone))/sum(area)-1/8)<2e-14, ...
    'The rhoC high-density region does not occupy one eighth of P16.');
rho(:,5)=coneScale*(1+(R-1)*double(cone));

for k=1:5
    for iz=1:nz
        jz=partner(iz);
        if jz>iz
            v=(rho(iz,k)+rho(jz,k))/2;
            rho([iz,jz],k)=v;
        end
    end
end

analyticPolygonMass=c0*template.polygonArea;
analyticFullDiskMass=4;
individualLower=[c0,rhoStar,rhoStar,rhoStar,coneScale];
contrast=[1,19,19,19,R];
cellRule={"analytic constant","analytic aligned step", ...
    "positive degree-4 triangle rule", ...
    "positive degree-4 triangle rule","analytic aligned step"};
parts=cell(1,5);
rows=struct([]);
for k=1:5
    p=template;
    p.rhoAvg=rho(:,k);
    p.mu=area.*p.rhoAvg;
    p.muLo=p.mu;
    p.muHi=p.mu;
    p.massSum=sum(p.mu);
    p.massRelRadius=0;
    p.massBoxes=ones(nz,1);
    p.kernel=char(names{k});
    p.kernelDescription=char(descriptions{k});
    p.beta=1/8;
    p.kernelMass=analyticFullDiskMass;
    p.kernelRstar=1;
    p.kernelRhoStar=individualLower(k);
    p.uniformRhoStar=rhoStar;
    p.kernelContrast=contrast(k);
    p.cellAverageRule=char(cellRule{k});
    p.fullDiskMassAnalytic=analyticFullDiskMass;
    p.fullDiskMassExactRule=full_disk_mass_rule(k,R,omega,phi);
    p.polygonMassAnalytic=analyticPolygonMass;
    if k==2, p.amplitude=0.9; end
    if k==3
        p.amplitudes=[0.45,0.45];
        p.phase2=pi/16;
        p.phase4=pi/32;
    elseif k==4
        p.amplitude=0.9;
        p.phase4=pi/32;
    elseif k==5
        p.contrastRatio=R;
        p.omega=omega;
        p.phi=phi;
        p.coneScale=coneScale;
    end

    evenRho=max(abs(p.rhoAvg-p.rhoAvg(partner)));
    evenMu=max(abs(p.mu-p.mu(partner)));
    polygonError=abs(p.massSum-analyticPolygonMass);
    diskError=abs(p.fullDiskMassExactRule-analyticFullDiskMass);
    assert(isequal(p.triangles,template.triangles));
    assert(p.qActual<=0.5*(1+32*eps),'The actual partition diameter exceeds 1/2.');
    assert(evenRho<5e-15 && evenMu<5e-14,'A kernel is not centrally even.');
    assert(all(p.rhoAvg>0) && min(p.rhoAvg)>=rhoStar-5e-14, ...
        'A kernel cell average violates the common positive lower bound.');
    assert(p.kernelRhoStar>=rhoStar-5e-14, ...
        'A pointwise kernel lower bound is below the common lower bound.');
    assert(polygonError<5e-13,'The five kernels do not have the same P16 mass.');
    assert(diskError<5e-13,'The exact-rule and analytic full-disk masses differ.');
    assert(p.kernelMass==4 && p.beta==1/8 && p.kernelRstar==1);
    parts{k}=p;

    row=struct( ...
        'kernel',names{k}, ...
        'description',descriptions{k}, ...
        'n_cells',nz, ...
        'q_actual',p.qActual, ...
        'cell_average_rule',cellRule{k}, ...
        'triangle_rule_degree',double(4*(k==3 || k==4)), ...
        'minimum_triangle_rule_weight',min(w), ...
        'reference_mass',p.kernelMass, ...
        'beta',p.beta, ...
        'r_star',p.kernelRstar, ...
        'kernel_rho_star',p.kernelRhoStar, ...
        'uniform_rho_star',p.uniformRhoStar, ...
        'contrast_bound',p.kernelContrast, ...
        'analytic_full_disk_mass',p.fullDiskMassAnalytic, ...
        'exact_rule_full_disk_mass',p.fullDiskMassExactRule, ...
        'full_disk_mass_error',diskError, ...
        'analytic_polygon_mass',p.polygonMassAnalytic, ...
        'computed_polygon_mass',p.massSum, ...
        'polygon_mass_error',polygonError, ...
        'cell_evenness_error',evenRho, ...
        'cell_mass_evenness_error',evenMu, ...
        'minimum_cell_average',min(p.rhoAvg), ...
        'cone_crossing_cells',double(k==5)*coneCrossings, ...
        'degree4_rule_moment_error',ruleError);
    if isempty(rows), rows=row; else, rows(end+1)=row; end %#ok<AGROW>
end
kernelQA=struct2table(rows);

polygonMasses=cellfun(@(p)p.massSum,parts);
assert(all(abs(polygonMasses-analyticPolygonMass)<5e-13));
for k=2:5
    assert(isequal(parts{1}.triangles,parts{k}.triangles), ...
        'All five kernels must share exactly the same triangles.');
end
end

function val=kernel_values(k,xy)
x=xy(:,1);
y=xy(:,2);
c0=4/pi;
switch k
    case 3
        a=pi/16;
        b=pi/32;
        h2=cos(2*a)*(x.^2-y.^2)+2*sin(2*a)*x.*y;
        h4=cos(4*b)*(x.^4-6*x.^2.*y.^2+y.^4)+ ...
            sin(4*b)*(4*x.*y.*(x.^2-y.^2));
        val=c0*(1+0.45*h2+0.45*h4);
    case 4
        b=pi/32;
        h4=cos(4*b)*(x.^4-6*x.^2.*y.^2+y.^4)+ ...
            sin(4*b)*(4*x.*y.*(x.^2-y.^2));
        val=c0*(1+0.9*h4);
    otherwise
        error('Polynomial values are implemented only for rhoG and rho4.');
end
assert(all(val>=0.4/pi-5e-14),'A polynomial kernel violates its lower bound.');
end

function [bary,w]=positive_degree4_rule()
a=0.816847572980459;
b=0.091576213509771;
c=0.108103018168070;
d=0.445948490915965;
bary=[a b b;b a b;b b a;c d d;d c d;d d c];
w=[repmat(0.109951743655322,3,1); ...
   repmat(0.223381589678011,3,1)];
end

function err=degree4_rule_error(bary,w)
x=bary(:,2);
y=bary(:,3);
err=0;
for i=0:4
    for j=0:4-i
        byRule=0.5*sum(w.*x.^i.*y.^j);
        exact=factorial(i)*factorial(j)/factorial(i+j+2);
        err=max(err,abs(byRule-exact));
    end
end
end

function mass=full_disk_mass_rule(k,R,omega,phi)
s=sqrt(3/5);
r=[(1-s)/2;1/2;(1+s)/2];
wr=[5/18;4/9;5/18];
nt=64;
theta=((0:nt-1)'+0.5)*(2*pi/nt);
x=cos(theta);
y=sin(theta);
c0=4/pi;
mass=0;
for ir=1:3
    xy=r(ir)*[x,y];
    switch k
        case 1
            value=c0*ones(nt,1);
        case 2
            value=c0*(1+0.9*sign(x.^2-y.^2));
        case {3,4}
            value=kernel_values(k,xy);
        case 5
            e=[cos(phi),sin(phi)];
            eperp=[-sin(phi),cos(phi)];
            inside=abs(xy*eperp')<tan(omega)*abs(xy*e');
            scale=4/(pi+2*omega*(R-1));
            value=scale*(1+(R-1)*double(inside));
        otherwise
            error('Unknown representative kernel.');
    end
    mass=mass+wr(ir)*r(ir)*(2*pi/nt)*sum(value);
end
end
