function [A,b] = cp_triangle_halfspace(v)
if size(v,1)~=3 || size(v,2)~=2
    error('v must be 3-by-2.');
end
area2 = det([v(2,:)-v(1,:);v(3,:)-v(1,:)]);
if area2<0, v=v([1 3 2],:); elseif area2==0, error('Degenerate triangle.'); end
A=zeros(3,2); b=zeros(3,1);
for k=1:3
    a=v(k,:); c=v(mod(k,3)+1,:); e=c-a;
    A(k,:)=[e(2),-e(1)];
    b(k,:)=A(k,:)*a';
end
end
