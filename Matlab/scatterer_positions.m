function [sp,g,w]=scatterer_positions(rotation,window,grid,center)
  %creates a Z^3 lattice of size grid/2, at center
  [x,y,z]=ndgrid((-grid/2:grid/2)+round(center(1)),(-grid/2:grid/2)+round(center(2)),(-grid/2:grid/2)+round(center(3)));
  g=[x(:),y(:),z(:)]';
  %rotates lattice
  p=rotation*g;
  %finds points within distance window about xy-plane, returns xy-coordinates
  w=find(abs(p(3,:))<window);
  sp=p(1:2,w);
end
