function out = mg_pcg_recomputed(A,b,Bfun,opts)

if nargin<4, opts = struct; end
if ~isfield(opts,'tol'), opts.tol = 1e-8; end
if ~isfield(opts,'maxit'), opts.maxit = 400; end

if isa(A,'function_handle')
    rawA = A;
else
    rawA = @(x) A*x;
end
if ~isa(Bfun,'function_handle')
    error('Bfun must apply an inverse-type preconditioner.');
end
b = b(:);
nA = 0;
nB = 0;
clock = tic;
[x,flag,reportedRelres,iter,resvec] = pcg(@countedA,b,opts.tol,opts.maxit,@countedB);
solveSeconds = toc(clock);
solverMatvecCalls = nA;
solverPreconditionerCalls = nB;
trueResidualVector = b-countedA(x); % deliberately independent recomputation
trueRelativeResidual = norm(trueResidualVector)/max(norm(b),realmin);

out = struct;
out.x = x;
out.flag = flag;
out.iterations = iter;
out.reportedRelativeResidual = reportedRelres;
out.recomputedRelativeResidual = trueRelativeResidual;
out.absoluteRecomputedResidual = norm(trueResidualVector);
out.solveSeconds = solveSeconds;
out.solverMatvecCalls = solverMatvecCalls;
out.recomputationMatvecCalls = 1;
out.totalMatvecCalls = nA;
out.preconditionerCalls = solverPreconditionerCalls;
out.residualHistory = resvec;
out.residualHistoryFinalRatio = resvec(end)/max(resvec(1),realmin);

    function y = countedA(xin)
        nA = nA+1;
        y = rawA(xin);
    end

    function y = countedB(xin)
        nB = nB+1;
        y = Bfun(xin);
    end
end
