function [sp,g,w]=scatterer_positions(rotation,window,pgrid,center,outdim,plot,radius)
  %creates a Z^dim lattice of size grid/2, at center
  dim = size(center,1);
  if dim == 3
    [x,y,z]=ndgrid((-pgrid/2:pgrid/2)+round(center(1)),(-pgrid/2:pgrid/2)+round(center(2)),(-pgrid/2:pgrid/2)+round(center(3)));
    g=[x(:),y(:),z(:)]';
  elseif dim == 4
    [x,y,z,h]=ndgrid((-pgrid/2:pgrid/2)+round(center(1)),(-pgrid/2:pgrid/2)+round(center(2)),(-pgrid/2:pgrid/2)+round(center(3)),(-pgrid/2:pgrid/2)+round(center(4)));
    g=[x(:),y(:),z(:),h(:)]';
  elseif dim == 5
    [x,y,z,h,l]=ndgrid((-pgrid/2:pgrid/2)+round(center(1)),(-pgrid/2:pgrid/2)+round(center(2)),(-pgrid/2:pgrid/2)+round(center(3)),(-pgrid/2:pgrid/2)+round(center(4)),(-pgrid/2:pgrid/2)+round(center(5)));
    g=[x(:),y(:),z(:),h(:),l(:)]';
  end
  p=rotation*g;
  w=find(abs(p(dim,:))<window);
  sp=p(1:outdim,w);
  if plot
    figure;
    hold on;
    numPoints = 10;
    % Loop through each center point
    for i = 1:size(sp, 2)
        center = sp(:,i);
        
        % Create a grid of points on the sphere's surface
        theta = linspace(0, 2*pi, numPoints);
        phi = linspace(0, pi, numPoints);
        [theta, phi] = meshgrid(theta, phi);
        
        x = radius * sin(phi) .* cos(theta) + center(1);
        y = radius * sin(phi) .* sin(theta) + center(2);
        z = radius * cos(phi) + center(3);
        
        % Plot the sphere's surface points
        surf(x, y, z, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    end
    
    % Set plot properties
    axis equal;
    grid on;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title('Spheres around Center Points');
    
    % Show the plot
    hold off;
  end
end
