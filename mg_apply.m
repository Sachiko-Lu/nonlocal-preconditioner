function x = mg_apply(mg,r,finestLevel)

if nargin<3 || isempty(finestLevel)
    finestLevel = numel(mg.level);
end
if size(r,1)~=mg.level{finestLevel}.nDof
    error('Right-hand side has the wrong number of rows.');
end
x = recurse(finestLevel,r);

    function z = recurse(ell,rhs)
        lev = mg.level{ell};
        if ell==1
            L = lev.coarseCholLower;
            z = L'\(L\rhs);
            return
        end

        z = zeros(size(rhs),'like',rhs);
        for sweep = 1:mg.nSmooth
            z = z+lev.omega*((rhs-lev.K*z)./lev.Dk);
        end

        residual = rhs-lev.K*z;
        P = mg.P{ell-1};
        z = z+P*recurse(ell-1,P'*residual);

        for sweep = 1:mg.nSmooth
            z = z+lev.omega*((rhs-lev.K*z)./lev.Dk);
        end
    end
end
