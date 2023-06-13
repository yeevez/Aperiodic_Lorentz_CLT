function paths=scatter_p(bounces,trials,step)
  m=[1,1,1;1,2,2;1,2,3];
  window=0.5;
  grid=20;
  ratio=0.99;
  [v,l]=eig(m);
  n=cross(v(:,2),v(:,3));
  %n=cross(v(:,1),v(:,3));
  u=cross(n,[0;0;1]);
  u=u/norm(u);
  t=acos(dot(n,[0;0;1]));
  r=[cos(t)+u(1)^2*(1-cos(t)),u(1)*u(2)*(1-cos(t))-u(3)*sin(t),u(1)*u(3)*(1-cos(t))+u(2)*sin(t);u(2)*u(1)*(1-cos(t))+u(3)*sin(t),cos(t)+u(2)^2*(1-cos(t)),u(2)*u(3)*(1-cos(t))-u(1)*sin(t);u(3)*u(1)*(1-cos(t))-u(2)*sin(t),u(3)*u(2)*(1-cos(t))+u(1)*sin(t),cos(t)+u(3)^2*(1-cos(t))];
  rinv=inv(r);
  [sp,~,~]=scatterer_positions(r,window,200,[0;0;0]);
  radius=ratio*min(arrayfun(@(i)min(sqrt(sum((sp(:,i)-sp(:,(i+1):end)).^2))),1:(size(sp,2)-1)))/2;
  fprintf('Using %f as the radius of the scatterers.\n',radius);
  paths=zeros(2,bounces/step+1,trials);
  parfor i=1:trials
    c=[0;0;0];
    [sp,~,~]=scatterer_positions(r,window,grid,c);
    angle=mod(2*pi*rand(),2*pi);
    position=radius*[cos(angle);sin(angle)];
    path=zeros(2,bounces+1);
    bounce=1;
    path(:,bounce)=position;
    tic;
    while bounce<=bounces
      rp=[cos(-angle),-sin(-angle);sin(-angle),cos(-angle)]*position;
      rsp=[cos(-angle),-sin(-angle);sin(-angle),cos(-angle)]*sp;
      h=find(rsp(2,:)>=(rp(2)-radius)&rsp(2,:)<=(rp(2)+radius)&rsp(1,:)>rp(1));
      if length(h)>0
        bounce=bounce+1;
        [~,s]=min(rsp(1,h));
        b=rsp(1,h(s))-sqrt(radius^2-(rp(2)-rsp(2,h(s)))^2);
        position=[cos(angle),-sin(angle);sin(angle),cos(angle)]*[b;rp(2)];
        angle=mod(2*atan2(rsp(2,h(s))-rp(2),rsp(1,h(s))-b)+pi+angle,2*pi);
        path(:,bounce)=position;
        %fprintf('Bounced at (%.2f,%.2f).\n',position(1),position(2));
      else
        h=rsp(:,convhull(rsp(1,:),rsp(2,:)));
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
        c=round(rinv*[[cos(angle),-sin(angle);sin(angle),cos(angle)]*[max(x(xi))-radius;rp(2)];0]);
        %fprintf('Regenerating grid at (%i,%i,%i).\n',c(1),c(2),c(3));
        [sp,~,~]=scatterer_positions(r,window,grid,c);
      end
    end
    e=toc;
    fprintf('Trial %i took %.2f seconds.\n',i,e);
    paths(:,:,i)=path(:,1:step:end);
  end
end
