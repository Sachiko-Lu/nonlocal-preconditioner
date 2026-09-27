function qa = mg_explicit_symmetric_spectrum(A,Bfun)

A = (A+A')/2;
n = size(A,1);
if size(A,2)~=n
    error('A must be square.');
end
Braw = full(Bfun(speye(n)));
symmetryError = norm(Braw-Braw','fro')/max(norm(Braw,'fro'),realmin);
B = (Braw+Braw')/2;
[L,p] = chol(B,'lower');
if p~=0
    error('The explicit inverse preconditioner is not positive definite.');
end
Sraw = L'*A*L;
SsymmetryError = norm(Sraw-Sraw','fro')/max(norm(Sraw,'fro'),realmin);
S = (Sraw+Sraw')/2;
[V,lambda] = eig(S,'vector');
[lambdaMin,imin] = min(real(lambda));
[lambdaMax,imax] = max(real(lambda));
if lambdaMin<=0
    error('The symmetric preconditioned matrix is not positive definite.');
end
resMin = norm(S*V(:,imin)-lambdaMin*V(:,imin))/max(abs(lambdaMin),realmin);
resMax = norm(S*V(:,imax)-lambdaMax*V(:,imax))/max(abs(lambdaMax),realmin);

qa = struct;
qa.nDof = n;
qa.preconditionerSymmetryError = symmetryError;
qa.minimumPreconditionerEigenvalue = min(eig(B));
qa.symmetricObjectSymmetryError = SsymmetryError;
qa.lambdaMin = lambdaMin;
qa.lambdaMax = lambdaMax;
qa.conditionNumber = lambdaMax/lambdaMin;
qa.extremeEigenResidual = max(resMin,resMax);
qa.spectralObject = 'L''*A*L, where B=L*L''';
end
