function simulation = simulation(matrix_list,bounces, trials, step, radius, n_existing_matrix_samples,outdim)
    %update this if we run again with more matrices, purely for file naming purposes
    for idx = 1:numel(matrix_list)
        matrix_number = idx + n_existing_matrix_samples;
        matrix = matrix_list{idx};
        [paths,max_fpl,scatterers] = scatter(bounces,trials,step,matrix,radius,outdim);
        if (length(paths) == 1)
            simulation = paths*0.99;
            return
        end
        save("scatter_samples_matrix" + matrix_number + ".mat" ,'paths');
        save("scatter_samples_matrix" +matrix_number + "_mfpl.mat","max_fpl");
        save("scatterer_positions_matrix" + matrix_number + ".mat" ,'scatterers');
    end
    simulation = -1;
end