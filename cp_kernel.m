function ker = cp_kernel(id)
id = upper(string(id));
switch id
    case "K1"
        ker.rho = @(r) (4/pi)*double(r<1 | abs(r-1)<32*eps);
        ker.mass = 4;
        ker.rStar = 1;
        ker.rhoStar = 4/pi;
    case "K2"
        ker.rho = @(r) (20/pi)*max(0,1-r);
        ker.mass = 20/3;
        ker.rStar = 1/2;
        ker.rhoStar = 10/pi;
    case "K3"
        ker.rho = @(r) (3/pi)./max(r,realmin);
        ker.mass = 6;
        ker.rStar = 1;
        ker.rhoStar = 3/pi;
    otherwise
        error('Unknown kernel %s.',id);
end
ker.id = char(id);
ker.beta = 1/(2*ker.mass);
end
