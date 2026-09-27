function spectrum = estimate_practical_spectrum_arnoldi(A,mg,delta,ell,opts)

if nargin < 5, opts = struct; end
if ~isfield(opts,'tol'), opts.tol = 3e-10; end
if ~isfield(opts,'maxit'), opts.maxit = 2000; end
if ~isfield(opts,'p'), opts.p = 64; end
if ~isfield(opts,'seed'), opts.seed = 20260821; end
if ~isfield(opts,'beta'), opts.beta = 1/8; end

n = size(A,1);
assert(size(A,2)==n && n==mg.level{ell}.nDof);
matrixSymmetryError = norm(A-A','fro')/max(norm(A,'fro'),realmin);
assert(matrixSymmetryError < 1e-11);
Bfun = @(v)mg_apply_practical(mg,v,opts.beta,delta,ell);
Tfun = @(v)Bfun(A*v);

oldState = rng;
cleanup = onCleanup(@()rng(oldState));
rng(opts.seed,'twister');
v0 = randn(n,1);
v0 = v0/norm(v0);
eigsOptions = struct('isreal',true,'issym',false,'tol',opts.tol, ...
    'maxit',opts.maxit,'p',min(n,max(20,opts.p)),'v0',v0);

clock = tic;
[vmin,dmin,flagMin] = eigs(Tfun,n,1,'smallestreal',eigsOptions);
eigsOptions.v0 = flipud(v0);
[vmax,dmax,flagMax] = eigs(Tfun,n,1,'largestreal',eigsOptions);
lambdaMin = real(dmin(1,1));
lambdaMax = real(dmax(1,1));
assert(flagMin==0 && flagMax==0,'The Arnoldi endpoint solve did not converge.');
assert(lambdaMin>0 && lambdaMax>=lambdaMin);

relativeResidualMin = norm(Tfun(vmin)-lambdaMin*vmin)/ ...
    max(abs(lambdaMin)*norm(vmin),realmin);
relativeResidualMax = norm(Tfun(vmax)-lambdaMax*vmax)/ ...
    max(abs(lambdaMax)*norm(vmax),realmin);
imaginaryPart = max(abs(imag([dmin(1,1),dmax(1,1)])));

u = v0;
v = circshift(v0,1);
Bu = Bfun(u);
Bv = Bfun(v);
preconditionerSymmetryError = abs(u'*Bv-v'*Bu)/ ...
    max(norm(u)*norm(Bv)+norm(v)*norm(Bu),realmin);

spectrum = struct( ...
    'lambdaMin',lambdaMin, ...
    'lambdaMax',lambdaMax, ...
    'conditionNumber',lambdaMax/lambdaMin, ...
    'minimumRelativeResidual',relativeResidualMin, ...
    'maximumRelativeResidual',relativeResidualMax, ...
    'extremeEigenResidual',max(relativeResidualMin,relativeResidualMax), ...
    'maximumImaginaryPart',imaginaryPart, ...
    'minimumFlag',flagMin, ...
    'maximumFlag',flagMax, ...
    'preconditionerSymmetryError',preconditionerSymmetryError, ...
    'matrixSymmetryError',matrixSymmetryError, ...
    'seconds',toc(clock), ...
    'method',"Arnoldi eigenvalues of B_h A_h");
clear cleanup
end
