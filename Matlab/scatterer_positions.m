function [sp,g,w]=scatterer_positions(rotation,window,grid,center)
  %creates a Z^dim lattice of size grid/2, at center
  dim = size(center);
  X=ndgrid_variable_dimensions(center,grid);
  g = cellfun(@(x) x(:), X)';
  p=rotation*g;
  w=find(abs(p(dim,:))<window);
  sp=p(1:outdim,w);
end
