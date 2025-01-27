function M = compute_M_parallel(masses, inertiae, jacobians)
%COMPUTE_M computes the generalized inertia matrix
%   The matrix M appears in the standard form of dynamics equation. It is
%   computed using jacobians of rigid bodies' CoM of the robot
    m = length(masses);             % # of rigid bodies
    n = length(jacobians(1,:)) / m; % # of configuration variables
    M = zeros(n);
    parfor i = 1:m
        J_i = jacobians(:,  n*(i-1)+1:n*i);
        I_i = inertiae(:, 3*(i-1)+1:3*i);
        m_i = masses(i);
        J_lin_i = J_i(1:3, :);
        J_rot_i = J_i(4:6, :);
        M = M + m_i * (J_lin_i' * J_lin_i) + (J_rot_i' * I_i * J_rot_i);
    end
end