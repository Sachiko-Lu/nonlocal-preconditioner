function z = mg_apply_practical(mg,r,beta,delta,finestLevel)
if nargin<5 || isempty(finestLevel)
    finestLevel = numel(mg.level);
end
if ~(isscalar(beta) && isfinite(beta) && beta>0)
    error('beta must be positive and finite.');
end
if ~(isscalar(delta) && isfinite(delta) && delta>0)
    error('delta must be positive and finite.');
end
D = mg.D{finestLevel};
z = mg_apply(mg,r,finestLevel)+(beta*delta^2)*(r./D);
end
