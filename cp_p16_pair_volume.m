function volume=cp_p16_pair_volume(delta)
theta=(0:15)'*pi/8;
P=delta*[cos(theta),sin(theta)];
P=clip(P,[1,0],1); P=clip(P,[-1,0],0);
P=clip(P,[0,1],1); P=clip(P,[0,-1],0);
lam=[2/3,1/6,1/6;1/6,2/3,1/6;1/6,1/6,2/3];
integral=0; c=mean(P,1);
for j=1:size(P,1)
    T=[c;P(j,:);P(mod(j,size(P,1))+1,:)];
    area=abs(det([T(2,:)-T(1,:);T(3,:)-T(1,:)]))/2;
    xy=lam*T;
    integral=integral+area/3*sum((1-xy(:,1)).*(1-xy(:,2)));
end
volume=2*delta^2*(8*sin(pi/8))-4*integral;
assert(volume>0);
end

function Q=clip(P,a,b)
Q=zeros(0,2);
for j=1:size(P,1)
    p=P(j,:);q=P(mod(j,size(P,1))+1,:);
    fp=a*p'-b;fq=a*q'-b;
    ip=fp<=0;iq=fq<=0;
    if ip,Q(end+1,:)=p;end %#ok<AGROW>
    if ip~=iq,Q(end+1,:)=p+(q-p)*(-fp/(fq-fp));end %#ok<AGROW>
end
end
