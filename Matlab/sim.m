matrices = {[1,1,1;1,1,2;1,2,3],[1,1,2;1,2,3;2,3,6],[1,2,1;2,5,4;1,4,6],[1,2,2;2,5,6;2,6,9],[1,2,2;2,5,5;2,5,6],[1,1,1;1,2,3;1,3,6]};
bounces = 1000000;
trials = 2000;
step = 1000;
radius = 0.2655*0.99;
n_existing_matrix_samples = 0; %update this if we run again with more matrices, purely for file naming purposes

for idx = 1:numel(matrix_list)
    matrix_number = idx + n_existing_matrix_samples;
    matrix = matrix_list{idx};
    [paths,max_fpl] = scatter(bounces,trials,step,matrix,radius);
    save("scatter_samples_matrix" + matrix_number + ".mat" ,'paths');
    save("scatter_samples_matrix" +matrix_number + "_mfpl.mat","max_fpl");
end
