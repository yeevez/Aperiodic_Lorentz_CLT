function [paths,flights]=scatter(bounces,trials,step,matrix,radius,outdim)
  indim = size(matrix,1);
  flights = zeros(trials);
  origin = zeros(size(matrix(1,:)))';
  window=0.5;
  %optimal grid size appears to be around 20x20x20
  pgrid=20;
  [v,~]=eig(matrix);
  %rotation matrix taking normal vector to z-axis
  r=rotation_matrix(v,outdim);
  %computes scatterer positions of a large grid
  [sp,~,~]=scatterer_positions(r,window,3,origin,outdim,false,radius);
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
  paths=zeros(outdim,bounces/step+1,trials);
  if outdim == 2
      parfor i=1:trials
        max_flight = 0;
        c=[0;0;0];
        [sp,~,~]=scatterer_positions(r,window,pgrid,c,outdim,false,radius);
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
          if ~isempty(h)
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
            [sp,~,~]=scatterer_positions(r,window,pgrid,c,outdim,false,radius);
          end
        end
        e=toc;
        fprintf('Trial %i took %.2f seconds.\n',i,e);
        flights(i) = max_flight
        paths(:,:,i)=path(:,1:step:end);
      end
  else 
      for i=1:trials
        max_flight = 0;
        c=origin;
        last_scatterer = [0;0;0];
        [sp,~,~]=scatterer_positions(r,window,pgrid,c,outdim,false,radius);
        %initial position and angle on surface of scatterer at origin
        theta = mod(2*pi*rand(),2*pi);
        phi = mod(pi*rand(),pi);
        position=radius*[sin(phi)*cos(theta);sin(phi)*sin(theta);cos(phi)];
        trajectory = position/norm(position);
        path=zeros(3,bounces+1);
        bounce=1;
        path(:,bounce)=position;
        tic;
        while bounce<=bounces
          %rotates so scattering direction is in +x-axis, possibly not at origin
          %trajectory = [sin(phi)*cos(theta);sin(phi)*sin(theta);cos(phi)]
          rot = rotation_matrix_alt(trajectory);
          rp = rot*position;
          %rotates scatterers, not considering the scatterer at the most recent collision
          rsp = rot*sp;
          %finds scatterers with center within radius of y-coordinate of
          %rotated position, towards the right, excluding the most recent
          %scatterer
          h=find(sqrt((rsp(2,:)-rp(2)).^2+(rsp(3,:)-rp(3,:)).^2)<=radius&rsp(1,:)>rp(1));
          if ~isempty(h)
            bounce=bounce+1;
            %reflect off of closest scatterer
            [~,s]=min(rsp(1,h));
            %position of reflection
            b=rsp(1,h(s))-sqrt(radius^2-(rp(2)-rsp(2,h(s)))^2-(rp(3)-rsp(3,h(s)))^2);
            %unrotated position of reflection
            position_new = [b;rp(2);rp(3)];
            position_new = rot\position_new;
            flight = norm(position_new-position);
            max_flight = max([flight,max_flight]);
            position = position_new;
            %angle of reflection
            last_scatterer = rot\rsp(:,h(s));
            lsp = last_scatterer-position;
            lsp = lsp/norm(lsp);
            trajectory = trajectory - 2*dot(lsp,trajectory)*lsp;
            %trajectory
            %theta = atan2(position(2)-last_scatterer(2),position(1)-last_scatterer(1))+pi;
            %phi= acos((position(3)-last_scatterer(3))/radius);
            path(:,bounce)=position;
            %fprintf('Bounced at (%.2f,%.2f).\n',position(1),position(2));
          else
            %passes through without reflection, computes convex hull of center of scatterers
            h=rsp(:,convhull(rsp(1,:),rsp(2,:),rsp(3,:)));
            n = size(h,2)/3;
            h = reshape(h,[size(h,1),n,3]);
            %computes intersection of scattering direction with convex hull, taking right-most intersection to be exit point
            intersections = [];
            for j=1:n
                points=h(:,j,:);
                p1 = points(:,:,1);
                p2 = points(:,:,2);
                p3 = points(:,:,3);
                normal = cross(points(:,:,2)-points(:,:,1),points(:,:,3)-points(:,:,1));
                %normal = normal/norm(normal);
                d = dot(normal,points(:,:,1));
                if normal(1) ~= 0
                    x = (d-(normal(2)*rp(2))-(normal(3)*rp(3)))/normal(1);
                    p = [x;rp(2);rp(3)];
                    if in_triangle(p,p1,p2,p3)
                        intersections = [intersections p];
                    end
                end
            end
            [~,idx] = sort(intersections(1,:),"descend");
            sorted_intersections = intersections(:,idx);
            intersection = sorted_intersections(:,1);
            %undo both rotations to be back in R^n
            if indim == 4
                c=round(r\[rot\intersection;0]);
            elseif indim == 5
                c=round(r\[rot\intersection;0;0]);
            end
            %computes new grid centered at exit point
            [sp,~,~]=scatterer_positions(r,window,pgrid,c,outdim,false,radius);
          end
        end
        e=toc;
        fprintf('Trial %i took %.2f seconds.\n',i,e);
        flights(i) = max_flight;
        paths(:,:,i)=path(:,1:step:end);
      end
  end
end
