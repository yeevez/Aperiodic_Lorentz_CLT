function [paths,flights]=scatterHypercube(bounces,trials,step,m,sideSize,lowerDim)
  flights = zeros(trials);
  origin = zeros(length(m),1);
  window=0.5;
  %optimal grid size appears to be around 20x20x20
  grid=25;
  [v,~]=eig(m);
  %normal vector to an eigenplane
  %rotation matrix taking normal vector to z-axis
  r=rotation_matrix(v,lowerDim);
  %computes scatterer positions of a large grid
  [sp,~,~]=scatterer_positions(r,window,200,origin,lowerDim);
  %disp(length(sp))
  %check if size will produce overlapping scatterers  
  %min(arrayfun(@(i)max(abs((sp(:,i)-sp(:,(i+1):end)).^2)),1:(size(sp,2)-1)))
  % jic = arrayfun(@(i) min(max(abs(sp(:,i)-sp(:,(i+1):end)))) , 1:(size(sp,2)-1));
  if sideSize > min( arrayfun(@(i) min(max(abs(sp(:,i)-sp(:,(i+1):end)))) , 1:(size(sp,2)-1)))/2
      disp("this size and scatterer configuration induces overlapping scatterers, try again with a smaller size")
      NewSize = min( arrayfun(@(i) min(max(abs(sp(:,i)-sp(:,(i+1):end)))) , 1:(size(sp,2)-1)))/2;
      disp("New Min size:")
      disp(NewSize)
      paths = NewSize;
      return 
  end
  fprintf('Using %f as the size of the scatterers.\n',sideSize);
  paths=zeros(2,bounces/step+1,trials);
  parfor i=1:trials
    max_flight = 0;
    [sp,~,~]=scatterer_positions(r,window,grid,origin, lowerDim);
    %initial position and angle on surface of scatterer at origin
    % norm(vel)
    angle=2*pi*rand();
    %Making it so that it aint inside a collider
    %position=sideSize*[cos(angle);sin(angle)];
    %+sp1(int16(rand()*length(sp1)))
    vel = [cos(angle);sin(angle)];
    %vel = vel/norm(vel);
    %not valid, inside the square.
    position=sideSize*vel/2;
    path=zeros(2,bounces+1);
    bounce=1;
    path(:,bounce)=position;
    tic
    while bounce<=bounces
      %rotates so scattering direction is in +x-axis, possibly not at origin 
      rp=[cos(-angle),-sin(-angle);sin(-angle),cos(-angle)]*position;
      %rotates scatterers
      rsp=[cos(-angle),-sin(-angle);sin(-angle),cos(-angle)]*sp;
      %finds scatterers with center within radius of y-coordinate of rotated position, towards the right
      radius = sideSize/sqrt(2);
      %finds scatterers with center within radius of y-coordinate of rotated position, towards the right
      %collision boolean and point of collision in side, position of reflection      
      big = find(rsp(2,:)>=(rp(2)-radius)&rsp(2,:)<=(rp(2)+radius)&rsp(1,:)>rp(1));
      small = find(rsp(2,:)>=(rp(2)-sideSize/2)&rsp(2,:)<=(rp(2)+sideSize/2)&rsp(1,:)>rp(1));
      boolColls = arrayfun(@(m) collWithHyperCube(position, sp(:, m), vel, sideSize), big);
      h = find(boolColls ~= 0);
      [~,pc,nv] = arrayfun(@(m) collWithHyperCube(position, sp(:,big(m)), vel, sideSize), h, "UniformOutput",false);
      pc = cell2mat(pc);
      nv = cell2mat(nv);
      if (length(small) > length(h))
          disp("Ehh")
      end
     if ~isempty(h)
        bounce=bounce+1;
        %reflect off of closest scatterer
        [~,s]= min( pc(1,:)-position(1) );
        % [~,s]= min(arrayfun(@(j) pc(1,j)-position(1), length(pc)));
        %position of reflection
        %unrotated position of reflection
        position_new = pc(:,s);
        vel = nv(:,s);
        flight = norm(position_new-position);
        max_flight = max(flight,max_flight);
        position = position_new;
        path(:,bounce)=position;
        angle = atan2(vel(2), vel(1));
        %fprintf('Bounced at (%.2f,%.2f).\n',position(1),position(2));0
      else
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
        c=round(r\[[cos(angle),-sin(angle);sin(angle),cos(angle)]*[max(x(xi))-radius;rp(2)];0]);
        [sp,~,~]=scatterer_positions(r,window,grid,c,lowerDim);
      end
    end
    e=toc;
    fprintf('Trial %i took %.2f seconds.\n',i,e);
    flights(i) = max_flight
    paths(:,:,i)=path(:,1:step:end);
  end
end