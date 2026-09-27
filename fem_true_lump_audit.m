function qa = fem_true_lump_audit(M,D,explicitLimit)
if nargin<3, explicitLimit = 2500; end
n = size(M,1);
Dlump = spdiags(D,0,n,n);
if n<=explicitLimit
    lambda = eig(full(Dlump),full(M));
    lambdaMin = min(real(lambda));
    lambdaMax = max(real(lambda));
    method = 'full generalized eigenvalue decomposition';
else
    opts.isreal = true; opts.tol = 1e-9;
    opts.maxit = 1000; opts.p = min(n,40);
    lambdaMin = eigs(Dlump,M,1,'smallestreal',opts);
    lambdaMax = eigs(Dlump,M,1,'largestreal',opts);
    method = 'symmetric generalized Lanczos';
end
qa = struct('nDof',n,'lambdaMinDrelativeM',real(lambdaMin), ...
    'lambdaMaxDrelativeM',real(lambdaMax), ...
    'lowerBoundViolation',max(0,1-real(lambdaMin)), ...
    'upperBoundViolation',max(0,real(lambdaMax)-4), ...
    'theoreticalInterval','[1,4] for 2D P1 true lumping','method',method);
end
