function [paths,flights]=scatter(bounces,trials,step,m,radius,outdim)
  flights = zeros(trials);
  origin = zeros(size(m(1,:)))';
  window=0.5;
  %optimal grid size appears to be around 20x20x20
  grid=20;
  [v,~]=eig(m);
  %rotation matrix taking normal vector to z-axis
  r=rotation_matrix(v,outdim);
  %rinv=inv(r);
  %computes scatterer positions of a large grid
  [sp,~,~]=scatterer_positions(r,window,20,origin,outdim);
  %check if radius will produce overlapping scatterers
  if radius > min(arrayfun(@(i)min(sqrt(sum((sp(:,i)-sp(:,(i+1):end)).^2))),1:(size(sp,2)-1)))/2
      disp("this radius and scatterer configuration induces overlapping scatterers, try again with a smaller radius")
      NewRadius = min(arrayfun(@(i)min(sqrt(sum((sp(:,i)-sp(:,(i+1):end)).^2))),1:(size(sp,2)-1)))/2;
      disp("New Min radius:")
      disp(NewRadius)
      paths = NewRadius;
      return 
  end
  fprintf('Using %f as the radius of the scatterers.\n',radius);
  paths=zeros(2,bounces/step+1,trials);
  if outdim == 2
      parfor i=1:trials
        max_flight = 0;
        c=[0;0;0];
        [sp,~,~]=scatterer_positions(r,window,grid,c,outdim);
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
            position_new =[cos(angle),-sin(angle);sin(angle),cos(angle)]*[b;rp(2)];
            flight = norm(position_new-position);
            max_flight = max([flight,max_flight])
            position = position_new
            %angle of reflection
            angle=mod(2*atan2(rsp(2,h(s))-rp(2),rsp(1,h(s))-b)+pi+angle,2*pi);
            path(:,bounce)=position;
            %fprintf('Bounced at (%.2f,%.2f).\n',position(1),position(2));
          else
            %passes through without reflection, computes convex hull of center of scatterers
            h=rsp(:,convhull(rsp(1,:),rsp(2,:)))
            %computes intersection of scattering direction with convex hull, taking right-most intersection to be exit point
            m=h(:,2:end)-h(:,1:(end-1));
            m=m(2,:)./m(1,:);
            x=(rp(2)-h(2,1:(end-1))+m.*h(1,1:(end-1)))./m;
            sh=sort([h(1,1:(end-1));h(1,2:end)],1);
            xi=find(x>sh(1,:)&x<sh(2,:));
            %undo both rotations to be back in R^3
            c=round(r\[[cos(angle),-sin(angle);sin(angle),cos(angle)]*[max(x(xi))-radius;rp(2)];0]);
            %fprintf('Regenerating grid at (%i,%i,%i).\n',c(1),c(2),c(3));
            %computes new grid centered at exit point
            [sp,~,~]=scatterer_positions(r,window,grid,c,outdim);
          end
        end
        e=toc;
        fprintf('Trial %i took %.2f seconds.\n',i,e);
        flights(i) = max_flight
        paths(:,:,i)=path(:,1:step:end);
      end
  else 
      parfor i=1:trials
        max_flight = 0;
        c=[0;0;0];
        [sp,~,~]=scatterer_positions(r,window,grid,c,outdim);
        %initial position and angle on surface of scatterer at origin
        theta =mod(2*pi*rand(),2*pi);
        phi = mod(pi*rand(),pi);
        position=radius*[sin(phi)*cos(theta);sin(phi)*sin(theta),cos(phi)];
        path=zeros(2,bounces+1);
        bounce=1;
        path(:,bounce)=position;
        tic;
        while bounce<=bounces
          %rotates so scattering direction is in +x-axis, possibly not at origin
          rot = rotation_matrix_alt(position);
          %%irot = inv(rot);
          rp = rot*position;
          %%rp=[cos(-angle),-sin(-angle);sin(-angle),cos(-angle)]*position;
          %rotates scatterers
          rsp = rot*sp
          %%rsp=[cos(-angle),-sin(-angle);sin(-angle),cos(-angle)]*sp;
          %finds scatterers with center within radius of y-coordinate of rotated position, towards the right
          h=find(norm([rsp(2,:);rsp(3,:)],2)<=radius&rsp(1,:)>rp(1));
          %%h=find(rsp(2,:)>=(rp(2)-radius)&rsp(2,:)<=(rp(2)+radius)&rsp(1,:)>rp(1));
          if length(h)>0
            bounce=bounce+1;
            %reflect off of closest scatterer
            [~,s]=min(rsp(1,h));
            %position of reflection
            b = rsp(1,h(s)) - sqrt(radius^2-(rp(2)-rsp(2,h(s)))^2-(rp(3)-rsp(3,h(s)))^2)
            %%b=rsp(1,h(s))-sqrt(radius^2-(rp(2)-rsp(2,h(s)))^2);
            %unrotated position of reflection
            position_new = [b;rp(2),rp(3)]\rot;
            %%position_new =[cos(angle),-sin(angle);sin(angle),cos(angle)]*[b;rp(2)];
            flight = norm(position_new-position);
            max_flight = max([flight,max_flight])
            position = position_new
            %angle of reflection
            theta=mod(2*atan2(rsp(2,h(s))-rp(2),rsp(1,h(s))-b)+pi+theta,2*pi);
            phi=mod(2*atan2(rsp(2,h(s))-rp(2),rsp(1,h(s))-b)+pi+phi,pi);
            %%angle=mod(2*atan2(rsp(2,h(s))-rp(2),rsp(1,h(s))-b)+pi+angle,2*pi);
            path(:,bounce)=position;
            %fprintf('Bounced at (%.2f,%.2f).\n',position(1),position(2));
          else
            %passes through without reflection, computes convex hull of center of scatterers
            psc = convhull(rsp(1,:),rsp(2,:),rsp(3,:));
            ps = []
            for j = 1:size(psc,1)
                ps = [ps;rsp(:,psc)]
            end
            

            h=rsp(:,convhull(rsp(1,:),rsp(2,:),rsp(3,:)));
            
            %computes intersection of scattering direction with convex hull, taking right-most intersection to be exit point
            m=h(:,2:end)-h(:,1:(end-1)); %todo: generalize % h = convex hull of scatterers -> m is columns 2:-1 of h - columns 1:-2 ????? what mean?
            m=m(2,:)./m(1,:); %todo: generalize %m is now the second row of m divided by first row -> y/x? slope?
            x=(rp(2)-h(2,1:(end-1))+m.*h(1,1:(end-1)))./m; %todo: generalize %pos_y -
            sh=sort([h(1,1:(end-1));h(1,2:end)],1); %todo: generalize
            xi=find(x>sh(1,:)&x<sh(2,:)); %todo: generalize
            %undo both rotations to be back in R^n
            c=round(r\[[cos(angle),-sin(angle);sin(angle),cos(angle)]*[max(x(xi))-radius;rp(2)];0]); %todo: generalize
            %computes new grid centered at exit point
            [sp,~,~]=scatterer_positions(r,window,grid,c,outdim);
          end
        end
        e=toc;
        fprintf('Trial %i took %.2f seconds.\n',i,e);
        flights(i) = max_flight
        paths(:,:,i)=path(:,1:step:end);
      end
  end
end
