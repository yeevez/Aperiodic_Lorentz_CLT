function R = rotation_matrix_alt(v)
    ax = zeros(size(v));
    for i = 1:size(v,2)
        ax(i,i) = 1;
    end
    C = v*ax';
    [U, ~, Vt] = svd(C);
    R = Vt * U';

    if R*v(:,1) < 0
    % If v goes to the negative x axis, rotate by an additional pi radians
    % so it goes to the positive x axis
    [vecs,vals] = eig(R);
    [~, idx] = min(abs(diag(vals) - 1));
    u = vecs(:,idx);
    u = u/norm(u);
    K = [0,-u(3),u(2);u(3),0,-u(1);-u(2),u(1),0];
    R_prime = eye(3) + 2*K^2; %rodrigues rotation formula 
    
    R = R_prime*R;
    end

end