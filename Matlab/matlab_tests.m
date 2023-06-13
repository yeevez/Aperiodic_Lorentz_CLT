c=[0;0;0];
m=[1,1,1;1,2,2;1,2,3];
[v,l]=eig(m);
n=cross(v(:,2),v(:,3));
u=cross(n,[0;0;1]);
u=u/norm(u);
t=acos(dot(n,[0;0;1]));
r=[cos(t)+u(1)^2*(1-cos(t)),u(1)*u(2)*(1-cos(t))-u(3)*sin(t),u(1)*u(3)*(1-cos(t))+u(2)*sin(t);u(2)*u(1)*(1-cos(t))+u(3)*sin(t),cos(t)+u(2)^2*(1-cos(t)),u(2)*u(3)*(1-cos(t))-u(1)*sin(t);u(3)*u(1)*(1-cos(t))-u(2)*sin(t),u(3)*u(2)*(1-cos(t))+u(1)*sin(t),cos(t)+u(3)^2*(1-cos(t))];
[sp,g,w] = scatterer_positions(r,2,6,c);

scatter(sp(1,:),sp(2,:))

save('matlabPoints.mat','sp')
