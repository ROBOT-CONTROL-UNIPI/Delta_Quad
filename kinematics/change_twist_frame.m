function L_p_q = change_twist_frame(Rpq)
% matrix for frame changing from {Q} to {P} maintaining the same pole

L_p_q = [     Rpq, zeros(3);
         zeros(3),      Rpq];

end