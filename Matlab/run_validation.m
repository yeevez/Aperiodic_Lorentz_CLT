load('scatter_samples_matrix1.mat', 'paths')
load('scatterer_positions_matrix1.mat', 'scatterers')
load('scatterer_generators1.mat', 'generators')
radius = .0850;
[vec,vec2] = validate(paths,scatterers,generators,radius);
vec
vec2