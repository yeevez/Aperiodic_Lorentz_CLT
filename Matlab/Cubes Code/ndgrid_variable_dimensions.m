function X = ndgrid_variable_dimensions(c, grid)
    num_dims = numel(c);
    ranges = cell(1, num_dims);
    for i = 1:num_dims
        center_i = c(i);
        range_i = (-grid/2:grid/2) + round(center_i);
        ranges{i} = range_i;
    end
    X = ndgrid(ranges{:});
end