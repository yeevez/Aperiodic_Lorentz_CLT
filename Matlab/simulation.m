function simulation = simulation(matrix_list,bounces, trials, step, radius, n_existing_matrix_samples)
    n_existing_matrix_samples = 0; 
    %update this if we run again with more matrices, purely for file naming purposes
    for idx = 1:numel(matrix_list)
        matrix_number = idx + n_existing_matrix_samples;
        matrix = matrix_list{idx};
        paths = scatter(bounces,trials,step,matrix,radius);
        if (length(paths) == 1)
            simulation = paths*0.99;
            return
        end
        save("scatter_samples_matrix" + matrix_number + ".mat" ,'paths');
    end
    simulation = -1;
end