load('scatter_samples_matrix1.mat', 'paths')
load('scatterer_positions_matrix1.mat', 'scatterers')

radius = .0850;
vec = validate(paths,scatterers,radius)