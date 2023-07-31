function X = ndgrid_variable_dimensions(c, grid)
    num_dims = numel(c);
    ranges = cell(1, num_dims);
    
    if num_dims == 1
        range = (-grid/2:grid/2) + round(c);
        ranges{1} = range;
    else
        for i = 1:num_dims
            center_i = c(i);
            range_i = (-grid/2:grid/2) + round(center_i);
            ranges{i} = range_i;
        end
        X = ndgrid(ranges{:});
    end
end