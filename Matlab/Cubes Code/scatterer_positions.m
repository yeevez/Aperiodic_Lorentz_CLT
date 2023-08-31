function [sp,g,w]=scatterer_positions(rotation,window,grid,center,outdim)
  %creates a Z^dim lattice of size grid/2, at center
  dim = size(center,1);
  if dim == 3
    [x,y,z]=ndgrid((-grid/2:grid/2)+round(center(1)),(-grid/2:grid/2)+round(center(2)),(-grid/2:grid/2)+round(center(3)));
    g=[x(:),y(:),z(:)]';
  elseif dim == 4
    [x,y,z,h]=ndgrid((-grid/2:grid/2)+round(center(1)),(-grid/2:grid/2)+round(center(2)),(-grid/2:grid/2)+round(center(3)),(-grid/2:grid/2)+round(center(4)));
    g=[x(:),y(:),z(:),h(:)]';
  elseif dim == 5
    [x,y,z,h,l]=ndgrid((-grid/2:grid/2)+round(center(1)),(-grid/2:grid/2)+round(center(2)),(-grid/2:grid/2)+round(center(3)),round(center(4)),round(center(5)));
    g=[x(:),y(:),z(:),h(:),l(:)]';
  end
  p=rotation*g;
  w=find(abs(p(dim,:))<window);
  sp=p(1:outdim,w);
end
