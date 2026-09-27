function mg = mg_setup(levels,Kcell,opts)

if nargin<2, Kcell = []; end
if nargin<3, opts = struct; end
opts = defaults(opts);
if ~(isscalar(opts.nSmooth) && opts.nSmooth>=1 && opts.nSmooth==round(opts.nSmooth))
    error('opts.nSmooth must be a positive integer.');
end
if ~(isscalar(opts.tau) && opts.tau>0 && opts.tau<2)
    error('opts.tau must lie strictly between 0 and 2.');
end

if isfield(levels,'mesh') && isfield(levels,'P')
    hierarchy = levels;
else
    error('levels must contain mesh and P cell arrays.');
end
nLevel = numel(hierarchy.mesh);
if numel(hierarchy.P)~=nLevel-1
    error('A hierarchy with L levels must contain L-1 prolongations.');
end

Kassembled = cell(nLevel,1);
M = cell(nLevel,1);
D = cell(nLevel,1);
massAudit = repmat(struct('true_lump_identity_error',NaN, ...
    'reduced_row_sum_relative_difference',NaN),nLevel,1);
for ell = 1:nLevel
    fem = fem_assemble_p1_operators(hierarchy.mesh{ell});
    Kassembled{ell} = fem.K;
    M{ell} = fem.M;
    D{ell} = fem.D;
    massAudit(ell).true_lump_identity_error = fem.trueLumpIdentityError;
    massAudit(ell).reduced_row_sum_relative_difference = ...
        fem.reducedRowSumRelativeDifference;
end

Kinput = normalize_Kcell(Kcell,Kassembled,nLevel);
K = Kinput;
if opts.useGalerkin
    for ell = nLevel-1:-1:1
        P = hierarchy.P{ell};
        K{ell} = sym(P'*K{ell+1}*P);
    end
end

level = cell(nLevel,1);
rediscretizationError = zeros(nLevel,1);
inputError = zeros(nLevel,1);
for ell = 1:nLevel
    A = sym(K{ell});
    n = size(A,1);
    if size(A,2)~=n || n~=numel(hierarchy.mesh{ell}.free)
        error('K{%d} has a dimension incompatible with the mesh.',ell);
    end
    dk = full(diag(A));
    if any(~isfinite(dk) | dk<=0)
        error('The Jacobi diagonal is not positive at level %d.',ell-1);
    end
    lambdaMax = scaled_largest_eigenvalue(A,dk,opts);
    omega = opts.tau/lambdaMax;
    level{ell} = struct('K',A,'Dk',dk,'M',M{ell},'D',D{ell}, ...
        'lambdaMaxDinvK',lambdaMax,'omega',omega,'nDof',n);
    rediscretizationError(ell) = relative_fro(A,Kassembled{ell});
    inputError(ell) = relative_fro(A,Kinput{ell});
end

L0 = chol(level{1}.K,'lower');
level{1}.coarseCholLower = L0;

mg = struct;
mg.level = level;
mg.K = cellfun(@(s)s.K,level,'UniformOutput',false);
mg.M = M;
mg.D = D;
mg.P = hierarchy.P;
mg.mesh = hierarchy.mesh;
mg.nSmooth = opts.nSmooth;
mg.tau = opts.tau;
mg.useGalerkin = logical(opts.useGalerkin);
mg.rediscretizationRelativeFro = rediscretizationError;
mg.inputRelativeFro = inputError;
mg.massAudit = massAudit;
mg.options = opts;
end

function opts = defaults(opts)
if ~isfield(opts,'nSmooth'), opts.nSmooth = 3; end
if ~isfield(opts,'tau'), opts.tau = 1.7; end
if ~isfield(opts,'useGalerkin'), opts.useGalerkin = true; end
if ~isfield(opts,'eigsTol'), opts.eigsTol = 1e-10; end
if ~isfield(opts,'eigsMaxit'), opts.eigsMaxit = 800; end
end

function K = normalize_Kcell(Kcell,Kassembled,nLevel)
if isempty(Kcell)
    K = Kassembled;
elseif iscell(Kcell)
    if numel(Kcell)~=nLevel
        error('Kcell must have one matrix per mesh level.');
    end
    K = Kcell(:);
else
    K = Kassembled;
    K{nLevel} = Kcell;
end
for ell = 1:nLevel
    K{ell} = sym(K{ell});
end
end

function lambdaMax = scaled_largest_eigenvalue(A,d,opts)
n = size(A,1);
if n<=80
    lambdaMax = max(eig(full(A),diag(d)));
else
    eigopts.isreal = true;
    eigopts.tol = opts.eigsTol;
    eigopts.maxit = opts.eigsMaxit;
    eigopts.p = min(n,max(24,2));
    eigopts.v0 = ones(n,1)/sqrt(n);
    [~,val,flag] = eigs(A,spdiags(d,0,n,n),1,'largestreal',eigopts);
    if flag~=0
        error('Largest scaled Jacobi eigenvalue did not converge.');
    end
    lambdaMax = real(val(1,1));
end
if ~(isfinite(lambdaMax) && lambdaMax>0)
    error('Invalid largest eigenvalue of D_K^{-1}K.');
end
end

function e = relative_fro(A,B)
den = norm(B,'fro');
if den==0, e = norm(A-B,'fro'); else, e = norm(A-B,'fro')/den; end
end

function A = sym(A)
A = (A+A')/2;
end
