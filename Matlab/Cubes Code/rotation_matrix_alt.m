function R = rotation_matrix_alt(v)
    ax = zeros(size(v));
    for i = 1:size(v,2)
        ax(i,i) = 1;
    end
    C = v*ax';
    [U, ~, Vt] = svd(C);
    R = Vt * U';

end