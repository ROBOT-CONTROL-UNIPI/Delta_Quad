function C = compute_C(M, q, q_dot)
%COMPUTE_C computes the Coriolis and centripetal forces matrix
%   The matrix C appears in the standard form of dynamics equation. It is
%   computed using Christoffel's symbols.
    n = length(q); 
    C = sym(zeros(n));
    gamma = sym(zeros(n, 1));
    for i = 1:n
        for j = 1:n
            for k = 1:n
                % compute gamma_ijk
                gamma(k) = 1/2 * ( diff( M(i,j), q(k) ) + diff( M(i,k), q(j) ) - ...
                                   diff( M(j,k), q(i) ) );
             
                C(i,j) = C(i,j) + gamma(k) * q_dot(k);
            end
        end
    end
end