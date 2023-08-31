function simulation = simulationHypercube(matrix_list,bounces, trials, step, sideSize, lowerDim, n_existing_matrix_samples)
    n_existing_matrix_samples = 0; 
    %update this if we run again with more matrices, purely for file naming purposes
    for idx = 1:numel(matrix_list)
        matrix_number = idx + n_existing_matrix_samples;
        matrix = matrix_list{idx};
        [paths,max_fpl] = scatterHypercube(bounces,trials,step,matrix,sideSize, lowerDim);
        % scatterHypercube(bounces,trials,step,m,size,lowerDim)
        if (length(paths) == 1)
            simulation = paths*0.99;
            return
        end
        save("scatter_samples_matrix" + matrix_number + ".mat" ,'paths');
        save("scatter_samples_matrix" + matrix_number + "_mfpl.mat","max_fpl");
    end
    simulation = -1;
end