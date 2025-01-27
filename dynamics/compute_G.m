function G = compute_G(masses, jacobians, g)
%COMPUTE_G computes gravity vector
%   The vector G appears in the standard form of dynamics equation. It is
%   computed as derivative of potential energy 
    m = length(masses);             % # of rigid bodies
    n = length(jacobians(1,:)) / m; % # of configuration variables
    G = sym(zeros(n, 1));
    for i = 1:n
        for j = 1:m
            J_vj = jacobians(:, n*(j-1)+1:n*j);
            G(i) = G(i) - masses(j)*g'*J_vj(1:3, i);
        end
    end
end