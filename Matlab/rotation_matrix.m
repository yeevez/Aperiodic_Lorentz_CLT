function R = rotation_matrix(v,outdim)
    n = size(v,1);
    v_new = zeros(n,outdim);
    for i = 1:outdim
        v_new(:,i) = v(:,end-i+1);
    end
    N = null(v_new'); %need way to decide sign of N
    ax = zeros(size(N));
    for i = 1:size(N,2)
        ax(end-size(N,2)+i,i) = 1;
    end
    C = N*ax';
    [U, ~, Vt] = svd(C);
    R = Vt * U';

end