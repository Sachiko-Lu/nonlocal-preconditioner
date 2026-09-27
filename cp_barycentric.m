function lam = cp_barycentric(T,x)
B=[T';ones(1,3)];
lam=(B\[x';ones(1,size(x,1))])';
end
