function M_p_q = change_twist_pole(P, Q)
% matrix for pole changing from Q to P in the same reference system

PQ = Q - P;
M_p_q = [  eye(3), hat(PQ);
         zeros(3),  eye(3)];

end