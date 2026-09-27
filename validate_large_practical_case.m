function result = validate_large_practical_case(A,hierarchy,mg,delta,ell,name)

mesh = hierarchy.mesh{ell};
n = numel(mesh.free);
beta = 1/8;
matrixSymmetry = norm(A-A','fro')/max(norm(A,'fro'),realmin);
assert(size(A,1)==n && size(A,2)==n && matrixSymmetry<1e-11);
Afun = @(x)A*x;
Bfun = @(x)mg_apply_practical(mg,x,beta,delta,ell);

K = mg.K{ell};
M = mg.M{ell};
LM = chol(M,'lower');
xy = mesh.p(mesh.free,:);
f = sin(pi*xy(:,1)).*sin(2*pi*xy(:,2)) ...
    +0.2*sin(3*pi*xy(:,1)).*sin(pi*xy(:,2));
rhs(1) = struct('name',"smooth_mass_load",'b',M*f);

oldState = rng;
cleanup = onCleanup(@()rng(oldState));
rng(20260821,'twister');
xi = 2*(rand(n,1)>=0.5)-1;
MinvXi = LM'\(LM\xi);
xi = xi/sqrt(xi'*MinvXi);
rhs(2) = struct('name',"normalized_Rademacher",'b',xi);

eigsOptions = struct('isreal',true,'tol',1e-10,'maxit',1000, ...
    'p',min(n,40));
[v,~,flag] = eigs(K,M,1,'largestreal',eigsOptions);
assert(flag==0,'The highest local generalized mode did not converge.');
v = real(v);
v = v/sqrt(v'*(M*v));
rhs(3) = struct('name',"highest_local_generalized_mode",'b',M*v);

u = randn(n,1);
v = randn(n,1);
Bu = Bfun(u);
Bv = Bfun(v);
preconditionerSymmetry = abs(u'*Bv-v'*Bu)/ ...
    max([abs(u'*Bv),abs(v'*Bu),realmin]);
minimumSampleRayleigh = min(u'*Bu/(u'*u),v'*Bv/(v'*v));
assert(preconditionerSymmetry<1e-11 && minimumSampleRayleigh>0);

operatorApplySeconds = median_apply_time(Afun,ones(n,1),3);
preconditionerApplySeconds = median_apply_time(Bfun,ones(n,1),3);
rows = cell(3,1);
for j = 1:3
    out = mg_pcg_recomputed(Afun,rhs(j).b,Bfun, ...
        struct('tol',1e-8,'maxit',300));
    row = common_row(mesh,delta,beta,name);
    row.preconditioner = "symmetric_MG_plus_true_lump";
    row.rhs = rhs(j).name;
    row.pcg_flag = out.flag;
    row.pcg_iterations = out.iterations;
    row.pcg_reported_relative_residual = out.reportedRelativeResidual;
    row.recomputed_relative_residual = out.recomputedRelativeResidual;
    row.solve_seconds = out.solveSeconds;
    row.solver_matvec_calls = out.solverMatvecCalls;
    row.residual_recomputation_matvec_calls = out.recomputationMatvecCalls;
    row.total_matvec_calls = out.totalMatvecCalls;
    row.preconditioner_calls = out.preconditionerCalls;
    row.preconditioner_apply_seconds = preconditionerApplySeconds;
    row.preconditioner_sample_symmetry_error = preconditionerSymmetry;
    row.preconditioner_minimum_sample_rayleigh = minimumSampleRayleigh;
    rows{j} = row;
end
pcgTable = struct2table(vertcat(rows{:}));
assert(all(pcgTable.pcg_flag==0));
assert(all(pcgTable.recomputed_relative_residual<1e-8));

estimated = estimate_practical_spectrum_arnoldi(A,mg,delta,ell);
srow = common_row(mesh,delta,beta,name);
srow.preconditioner = "symmetric_MG_plus_true_lump";
srow.lambda_min = estimated.lambdaMin;
srow.lambda_max = estimated.lambdaMax;
srow.condition_number = estimated.conditionNumber;
srow.preconditioner_symmetry_error = estimated.preconditionerSymmetryError;
srow.preconditioner_min_eigenvalue = NaN;
srow.symmetric_object_symmetry_error = NaN;
srow.extreme_eigen_residual = estimated.extremeEigenResidual;
srow.spectral_object = estimated.method;
spectrumTable = struct2table(srow);

qa = struct('nDof',n,'level',ell-1,'delta',delta,'beta',beta, ...
    'operatorApplySeconds',operatorApplySeconds, ...
    'practicalPreconditionerApplySeconds',preconditionerApplySeconds, ...
    'assemblySymmetryErrorBeforeAnySymmetrization',matrixSymmetry, ...
    'spectral',estimated);
result = struct('pcg',pcgTable,'spectrum',spectrumTable,'qa',qa);
clear cleanup
end

function row = common_row(mesh,delta,beta,name)
row = struct;
row.case_name = string(name);
row.kernel = "K1";
row.n_dof = numel(mesh.free);
row.n_elements = size(mesh.t,1);
row.hmax = mesh.hmax;
row.hmin = mesh.hmin;
row.delta = delta;
row.delta_over_hmax = delta/mesh.hmax;
row.beta = beta;
row.m_polygon = 16;
row.q_kernel_cell = 0.5;
row.eta_Q = NaN;
row.operator_field = "Adbl";
row.geometry_certified = false;
row.mass_enclosure_is_formal = false;
end

function seconds = median_apply_time(fun,x,nRepeat)
sample = zeros(nRepeat,1);
fun(x);
for j = 1:nRepeat
    clock = tic;
    fun(x);
    sample(j) = toc(clock);
end
seconds = median(sample);
end
