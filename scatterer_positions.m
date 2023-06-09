function [sp,g,w]=scatterer_positions(rotation,window,grid,center)
  [x,y,z]=ndgrid((-grid/2:grid/2)+round(center(1)),(-grid/2:grid/2)+round(center(2)),(-grid/2:grid/2)+round(center(3)));
  g=[x(:),y(:),z(:)]';
  p=rotation*g;
  w=find(abs(p(3,:))<window);
  sp=p(1:2,w);
end