m = [1,1,2;1,2,3;2,3,6];

[v,l]=eig(m);
%normal vector to an eigenplane
% n=cross(v(:,2),v(:,3));
n=cross(v(:,1),v(:,3));
n=n/norm(n);
u=cross(n,[0;0;1]);
u=u/norm(u);
t=acos(dot(n,[0;0;1]));
%rotation matrix taking normal vector to z-axis
r=[cos(t)+u(1)^2*(1-cos(t)),u(1)*u(2)*(1-cos(t))-u(3)*sin(t),u(1)*u(3)*(1-cos(t))+u(2)*sin(t);u(2)*u(1)*(1-cos(t))+u(3)*sin(t),cos(t)+u(2)^2*(1-cos(t)),u(2)*u(3)*(1-cos(t))-u(1)*sin(t);u(3)*u(1)*(1-cos(t))-u(2)*sin(t),u(3)*u(2)*(1-cos(t))+u(1)*sin(t),cos(t)+u(3)^2*(1-cos(t))];
rinv=inv(r);
window = 0.5;
grid = 6;
center = [0;0;0];
%creates a Z^3 lattice of size grid/2, at center
[x,y,z]=ndgrid((-grid/2:grid/2)+round(center(1)),(-grid/2:grid/2)+round(center(2)),(-grid/2:grid/2)+round(center(3)));
g=[x(:),y(:),z(:)]';
%rotates lattice
fprintf('mdsamdas')
p=r*g
rotation = inv(v);
fprintf("smkldamkas")
lp = rotation*g
mp = v\g