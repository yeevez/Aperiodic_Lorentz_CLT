function paths=scatter(bounces,trials,step)
  m=[1,1,1;1,2,2;1,2,3];
  %non-pisot
  %m=[3,3,5;1,4,3;2,1,3];
  window=0.5;
  %optimal grid size appears to be around 20x20x20
  grid=20;
  ratio=0.99;
  [v,l]=eig(m);
  %normal vector to an eigenplane
  n=cross(v(:,2),v(:,3));
  %n=cross(v(:,1),v(:,3));
  n=n/norm(n);
  u=cross(n,[0;0;1]);
  u=u/norm(u);
  t=acos(dot(n,[0;0;1]));
  %rotation matrix taking normal vector to z-axis
  r=[cos(t)+u(1)^2*(1-cos(t)),u(1)*u(2)*(1-cos(t))-u(3)*sin(t),u(1)*u(3)*(1-cos(t))+u(2)*sin(t);u(2)*u(1)*(1-cos(t))+u(3)*sin(t),cos(t)+u(2)^2*(1-cos(t)),u(2)*u(3)*(1-cos(t))-u(1)*sin(t);u(3)*u(1)*(1-cos(t))-u(2)*sin(t),u(3)*u(2)*(1-cos(t))+u(1)*sin(t),cos(t)+u(3)^2*(1-cos(t))];
  rinv=inv(r);
  %computes scatterer positions of a large grid
  [sp,~,~]=scatterer_positions(r,window,200,[0;0;0]);
  %set radius to be ratio of minimum pairwise distance between scatterers
  radius=ratio*min(arrayfun(@(i)min(sqrt(sum((sp(:,i)-sp(:,(i+1):end)).^2))),1:(size(sp,2)-1)))/2;
  fprintf('Using %f as the radius of the scatterers.\n',radius);
  paths=zeros(2,bounces/step+1,trials);
  parfor i=1:trials
    c=[0;0;0];
    [sp,~,~]=scatterer_positions(r,window,grid,c);
    %initial position and angle on surface of scatterer at origin
    angle=mod(2*pi*rand(),2*pi);
    position=radius*[cos(angle);sin(angle)];
    path=zeros(2,bounces+1);
    bounce=1;
    path(:,bounce)=position;
    tic;
    while bounce<=bounces
      %rotates so scattering direction is in +x-axis, possibly not at origin
      rp=[cos(-angle),-sin(-angle);sin(-angle),cos(-angle)]*position;
      %rotates scatterers
      rsp=[cos(-angle),-sin(-angle);sin(-angle),cos(-angle)]*sp;
      %finds scatterers with center within radius of y-coordinate of rotated position, towards the right
      h=find(rsp(2,:)>=(rp(2)-radius)&rsp(2,:)<=(rp(2)+radius)&rsp(1,:)>rp(1));
      if length(h)>0
        bounce=bounce+1;
        %reflect off of closest scatterer
        [~,s]=min(rsp(1,h));
        %position of reflection
        b=rsp(1,h(s))-sqrt(radius^2-(rp(2)-rsp(2,h(s)))^2);
        %unrotated position of reflection
        position=[cos(angle),-sin(angle);sin(angle),cos(angle)]*[b;rp(2)];
        %angle of reflection
        angle=mod(2*atan2(rsp(2,h(s))-rp(2),rsp(1,h(s))-b)+pi+angle,2*pi);
        path(:,bounce)=position;
        %fprintf('Bounced at (%.2f,%.2f).\n',position(1),position(2));
      else
        %passes through without reflection, computes convex hull of center of scatterers
        h=rsp(:,convhull(rsp(1,:),rsp(2,:)));
        %computes intersection of scattering direction with convex hull, taking right-most intersection to be exit point
        m=h(:,2:end)-h(:,1:(end-1));
        m=m(2,:)./m(1,:);
        x=(rp(2)-h(2,1:(end-1))+m.*h(1,1:(end-1)))./m;
        sh=sort([h(1,1:(end-1));h(1,2:end)],1);
        %xi=find(x>rp(1)&x>sh(1,:)&x<sh(2,:));
        %if length(xi)>0
        %  position=[cos(angle),-sin(angle);sin(angle),cos(angle)]*[max(x(xi))-radius;rp(2)];
        %end
        %c=round(rinv*[position;0]);
        xi=find(x>sh(1,:)&x<sh(2,:));
        %undo both rotations to be back in R^3
        c=round(rinv*[[cos(angle),-sin(angle);sin(angle),cos(angle)]*[max(x(xi))-radius;rp(2)];0]);
        %fprintf('Regenerating grid at (%i,%i,%i).\n',c(1),c(2),c(3));
        %computes new grid centered at exit point
        [sp,~,~]=scatterer_positions(r,window,grid,c);
      end
    end
    e=toc;
    fprintf('Trial %i took %.2f seconds.\n',i,e);
    paths(:,:,i)=path(:,1:step:end);
  end
end
