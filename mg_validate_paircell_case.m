function result = mg_validate_paircell_case(nl,hierarchy,mg,opts)

if nargin<4, opts = struct; end
opts = defaults(opts,mg);
ell = opts.level;
n = mg.level{ell}.nDof;
[Aexplicit,Afun,matrixField,assemblySymmetry] = get_operator(nl,n,opts);
[delta,beta] = get_parameters(nl,opts);
alpha = beta*delta^2;
K = mg.K{ell}; M = mg.M{ell}; D = mg.D{ell};
mesh = hierarchy.mesh{ell};

if ~isempty(Aexplicit)
    [~,spdFlag] = chol(Aexplicit);
    if spdFlag~=0
        error('The explicit pair-cell matrix is not positive definite.');
    end
end

LK = [];
LM = chol(M,'lower');
if n<=opts.directLimit
    LK = chol(K,'lower');
end
solveM = @(r) LM'\(LM\r);
if ~isempty(LK)
    solveK = @(r) LK'\(LK\r);
end

names = strings(0,1);
bfuns = cell(0,1);
if ~opts.practicalOnly && ~isempty(LK)
    names(end+1,1) = "K_inverse";
    bfuns{end+1,1} = solveK;
    names(end+1,1) = "ideal_K_plus_consistent_mass";
    bfuns{end+1,1} = @(r)solveK(r)+alpha*solveM(r);
end
names(end+1,1) = "symmetric_MG_plus_true_lump";
bfuns{end+1,1} = @(r)mg_apply_practical(mg,r,beta,delta,ell);

rhs = build_rhs(mesh,K,M,LM,opts);
pcgRows = cell(0,1);
row = 0;
for ip = 1:numel(bfuns)
    Bfun = bfuns{ip};
    [symmetryError,minSample] = sampled_preconditioner_qa(Bfun,n,opts.seed+31*ip);
    preApplySeconds = median_apply_time(Bfun,ones(n,1),opts.timingRepeats);
    for ir = 1:numel(rhs)
        row = row+1;
        out = mg_pcg_recomputed(Afun,rhs(ir).b,Bfun, ...
            struct('tol',opts.tol,'maxit',opts.maxit));
        qrow = common_row(nl,mesh,n,delta,beta,matrixField);
        qrow.preconditioner = names(ip);
        qrow.rhs = rhs(ir).name;
        qrow.pcg_flag = out.flag;
        qrow.pcg_iterations = out.iterations;
        qrow.pcg_reported_relative_residual = out.reportedRelativeResidual;
        qrow.recomputed_relative_residual = out.recomputedRelativeResidual;
        qrow.solve_seconds = out.solveSeconds;
        qrow.solver_matvec_calls = out.solverMatvecCalls;
        qrow.residual_recomputation_matvec_calls = out.recomputationMatvecCalls;
        qrow.total_matvec_calls = out.totalMatvecCalls;
        qrow.preconditioner_calls = out.preconditionerCalls;
        qrow.preconditioner_apply_seconds = preApplySeconds;
        qrow.preconditioner_sample_symmetry_error = symmetryError;
        qrow.preconditioner_minimum_sample_rayleigh = minSample;
        pcgRows{row,1} = qrow;
    end
end
pcgTable = struct2table(vertcat(pcgRows{:}));

spectralRows = cell(0,1);
if ~isempty(Aexplicit) && n<=opts.spectralLimit
    for ip = 1:numel(bfuns)
        qaSpec = mg_explicit_symmetric_spectrum(Aexplicit,bfuns{ip});
        srow = common_row(nl,mesh,n,delta,beta,matrixField);
        srow.preconditioner = names(ip);
        srow.lambda_min = qaSpec.lambdaMin;
        srow.lambda_max = qaSpec.lambdaMax;
        srow.condition_number = qaSpec.conditionNumber;
        srow.preconditioner_symmetry_error = qaSpec.preconditionerSymmetryError;
        srow.preconditioner_min_eigenvalue = qaSpec.minimumPreconditionerEigenvalue;
        srow.symmetric_object_symmetry_error = qaSpec.symmetricObjectSymmetryError;
        srow.extreme_eigen_residual = qaSpec.extremeEigenResidual;
        srow.spectral_object = string(qaSpec.spectralObject);
        spectralRows{ip,1} = srow;
    end
    spectralTable = struct2table(vertcat(spectralRows{:}));
else
    spectralTable = table;
end

massQa = fem_true_lump_audit(M,D,opts.massAuditExplicitLimit);
u = sin((1:n)'*sqrt(2));
AapplySeconds = median_apply_time(Afun,u,opts.timingRepeats);
BpracApplySeconds = median_apply_time(bfuns{end},u,opts.timingRepeats);
explicitApplySeconds = NaN;
if ~isempty(Aexplicit)
    explicitApplySeconds = median_apply_time(@(x)Aexplicit*x,u,opts.timingRepeats);
end
matrixFreeApplySeconds = NaN;
if isfield(nl,'applyA') && isa(nl.applyA,'function_handle')
    matrixFreeApplySeconds = median_apply_time(nl.applyA,u,opts.timingRepeats);
end
interface = optional_interface_errors(nl,K,M,D);
matrixFreeAudit = explicit_vs_apply_audit(nl,Aexplicit,n,opts.seed);

qa = struct;
qa.nDof = n;
qa.level = ell-1;
qa.delta = delta;
qa.beta = beta;
qa.alpha = alpha;
qa.matrixField = matrixField;
qa.geometryCertified = logical(get_optional(nl,'geometryCertified',false));
qa.massEnclosureIsFormal = logical(get_optional(nl,'massEnclosureIsFormal',false));
qa.explicitMatrixAvailable = ~isempty(Aexplicit);
qa.assemblySymmetryErrorBeforeAnySymmetrization = assemblySymmetry;
qa.operatorApplySeconds = AapplySeconds;
qa.explicitSparseMatvecSeconds = explicitApplySeconds;
qa.matrixFreeBlockApplySeconds = matrixFreeApplySeconds;
qa.practicalPreconditionerApplySeconds = BpracApplySeconds;
qa.mass = massQa;
qa.interface = interface;
qa.matrixFree = matrixFreeAudit;
qa.recomputedResidualDefinition = 'norm(b-A*x,2)/norm(b,2), with one fresh A application after PCG';
qa.spectralObject = 'L''*A*L for B=L*L''; never singular values of B*A';
result = struct('pcg',pcgTable,'spectrum',spectralTable,'qa',qa);
end

function opts = defaults(opts,mg)
if ~isfield(opts,'level'), opts.level = numel(mg.level); end
if ~isfield(opts,'delta'), opts.delta = []; end
if ~isfield(opts,'beta'), opts.beta = []; end
if ~isfield(opts,'tol'), opts.tol = 1e-8; end
if ~isfield(opts,'maxit'), opts.maxit = 400; end
if ~isfield(opts,'spectralLimit'), opts.spectralLimit = 1400; end
if ~isfield(opts,'directLimit'), opts.directLimit = 6000; end
if ~isfield(opts,'highModeLimit'), opts.highModeLimit = 2500; end
if ~isfield(opts,'massAuditExplicitLimit'), opts.massAuditExplicitLimit = 2500; end
if ~isfield(opts,'seed'), opts.seed = 20260821; end
if ~isfield(opts,'timingRepeats'), opts.timingRepeats = 3; end
if ~isfield(opts,'maxAssemblySymmetryError'), opts.maxAssemblySymmetryError = 1e-10; end
if ~isfield(opts,'practicalOnly'), opts.practicalOnly = false; end
end

function [A,Afun,fieldName,symmetryError] = get_operator(nl,n,opts)
A = [];
fieldName = "matrix_free_applyA";
candidateFields = {'Adbl','Ahat','Aavg','A'};
for j = 1:numel(candidateFields)
    f = candidateFields{j};
    if isfield(nl,f) && ~isempty(nl.(f))
        Araw = nl.(f);
        if ~ismatrix(Araw) || any(size(Araw)~=[n n])
            error('nl.%s has the wrong size.',f);
        end
        symmetryError = norm(Araw-Araw','fro')/max(norm(Araw,'fro'),realmin);
        if symmetryError>opts.maxAssemblySymmetryError
            error('nl.%s is not symmetric before postprocessing (relative error %.3e).',f,symmetryError);
        end
        A = Araw;
        Afun = @(x)A*x;
        fieldName = string(f);
        return
    end
end
if isfield(nl,'applyA') && isa(nl.applyA,'function_handle')
    Afun = nl.applyA;
    symmetryError = sampled_operator_symmetry(Afun,n);
else
    error('nl must provide Adbl, Ahat, Aavg, A, or applyA.');
end
end

function e = sampled_operator_symmetry(Afun,n)
old = rng; cleanup = onCleanup(@()rng(old)); rng(8317,'twister');
u = randn(n,1); v = randn(n,1);
lhs = u'*Afun(v); rhs = v'*Afun(u);
e = abs(lhs-rhs)/max([abs(lhs),abs(rhs),realmin]);
clear cleanup
end

function [delta,beta] = get_parameters(nl,opts)
if ~isempty(opts.delta), delta = opts.delta;
elseif isfield(nl,'delta'), delta = nl.delta;
elseif isfield(nl,'stats') && isfield(nl.stats,'delta'), delta = nl.stats.delta;
else, error('delta is missing from both nl and opts.');
end
if ~isempty(opts.beta), beta = opts.beta;
elseif isfield(nl,'beta'), beta = nl.beta;
elseif isfield(nl,'betaStar'), beta = nl.betaStar;
elseif isfield(nl,'stats') && isfield(nl.stats,'beta'), beta = nl.stats.beta;
else, error('beta is missing from both nl and opts.');
end
if ~(isscalar(delta)&&isfinite(delta)&&delta>0&&isscalar(beta)&&isfinite(beta)&&beta>0)
    error('delta and beta must be positive finite scalars.');
end
end

function rhs = build_rhs(mesh,K,M,LM,opts)
free = mesh.free(:);
if isfield(mesh,'p'), xy = mesh.p(free,:); else, xy = mesh.nodes(free,:); end
n = numel(free);
f = sin(pi*xy(:,1)).*sin(2*pi*xy(:,2))+0.2*sin(3*pi*xy(:,1)).*sin(pi*xy(:,2));
rhs(1) = struct('name',"smooth_mass_load",'b',M*f);

old = rng; cleanup = onCleanup(@()rng(old)); rng(opts.seed,'twister');
xi = 2*(rand(n,1)>=0.5)-1;
MinvXi = LM'\(LM\xi);
xi = xi/sqrt(xi'*MinvXi);
rhs(2) = struct('name',"normalized_Rademacher",'b',xi);

if n<=opts.highModeLimit
    eopts.isreal = true; eopts.tol = 1e-10;
    eopts.maxit = 1000; eopts.p = min(n,40);
    [v,~,flag] = eigs(K,M,1,'largestreal',eopts);
    if flag~=0, error('Highest local generalized mode did not converge.'); end
    v = real(v); v = v/sqrt(v'*(M*v));
    rhs(3) = struct('name',"highest_local_generalized_mode",'b',M*v);
end
clear cleanup
end

function [e,minRayleigh] = sampled_preconditioner_qa(Bfun,n,seed)
old = rng; cleanup = onCleanup(@()rng(old)); rng(seed,'twister');
u = randn(n,1); v = randn(n,1);
Bu = Bfun(u); Bv = Bfun(v);
lhs = u'*Bv; rhs = v'*Bu;
e = abs(lhs-rhs)/max([abs(lhs),abs(rhs),realmin]);
minRayleigh = min(u'*Bu/(u'*u),v'*Bv/(v'*v));
if minRayleigh<=0, error('A sampled preconditioner Rayleigh quotient is nonpositive.'); end
clear cleanup
end

function seconds = median_apply_time(fun,x,nRepeat)
sample = zeros(nRepeat,1);
fun(x); % warm-up
for j = 1:nRepeat
    c = tic; fun(x); sample(j) = toc(c);
end
seconds = median(sample);
end

function row = common_row(nl,mesh,n,delta,beta,matrixField)
row = struct;
row.case_name = string(get_optional(nl,'name','unnamed_paircell_case'));
row.kernel = string(get_optional(nl,'kernel','unspecified'));
row.n_dof = n;
if isfield(mesh,'nElem'), row.n_elements = mesh.nElem; else, row.n_elements = size(mesh.t,1); end
row.hmax = mesh.hmax;
row.hmin = mesh.hmin;
row.delta = delta;
row.delta_over_hmax = delta/mesh.hmax;
row.beta = beta;
row.m_polygon = get_numeric(nl,'m',NaN);
row.q_kernel_cell = get_numeric(nl,'q',NaN);
row.eta_Q = get_numeric(nl,'etaQ',NaN);
row.operator_field = string(matrixField);
row.geometry_certified = logical(get_optional(nl,'geometryCertified',false));
row.mass_enclosure_is_formal = logical(get_optional(nl,'massEnclosureIsFormal',false));
end

function out = get_optional(s,f,default)
if isfield(s,f), out = s.(f); else, out = default; end
end

function out = get_numeric(s,f,default)
if isfield(s,f) && isscalar(s.(f)), out = double(s.(f));
elseif isfield(s,'stats') && isfield(s.stats,f) && isscalar(s.stats.(f)), out = double(s.stats.(f));
elseif strcmp(f,'etaQ') && isfield(s,'etaQReported') && isscalar(s.etaQReported)
    out = double(s.etaQReported);
else, out = default;
end
end

function qa = optional_interface_errors(nl,K,M,D)
qa = struct('KrelativeFro',NaN,'MrelativeFro',NaN,'Drelative2',NaN);
if isfield(nl,'K') && ~isempty(nl.K)
    qa.KrelativeFro = norm(nl.K-K,'fro')/max(norm(K,'fro'),realmin);
end
if isfield(nl,'M') && ~isempty(nl.M)
    qa.MrelativeFro = norm(nl.M-M,'fro')/max(norm(M,'fro'),realmin);
end
if isfield(nl,'D') && ~isempty(nl.D)
    qa.Drelative2 = norm(nl.D(:)-D)/max(norm(D),realmin);
elseif isfield(nl,'Mlump') && ~isempty(nl.Mlump)
    qa.Drelative2 = norm(nl.Mlump(:)-D)/max(norm(D),realmin);
end
end

function qa = explicit_vs_apply_audit(nl,A,n,seed)
qa = struct('available',false,'relativeError',NaN,'mode','unavailable');
if isempty(A) || ~isfield(nl,'applyA') || ~isa(nl.applyA,'function_handle')
    return
end
qa.available = true;
if n<=200
    X = speye(n);
    qa.relativeError = norm(nl.applyA(X)-A*X,'fro')/max(norm(A,'fro'),realmin);
    qa.mode = 'all coordinate vectors';
else
    old = rng; cleanup = onCleanup(@()rng(old)); rng(seed+9029,'twister');
    X = randn(n,3);
    qa.relativeError = norm(nl.applyA(X)-A*X,'fro')/max(norm(A*X,'fro'),realmin);
    qa.mode = 'three deterministic random vectors';
    clear cleanup
end
end
