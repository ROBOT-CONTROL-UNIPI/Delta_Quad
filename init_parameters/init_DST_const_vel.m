syms x_0 x_1 ts vel real
v_0 = -vel;
a_0 = 0;
v_1 = -vel;
a_1 = 0;

syms x_0 x1 x2 x3 x4 x5 x6 t real
p = x1 + x2*t + x3*t^2 + x4*t^3 +x5*t^4 + x6*t^5;
p_dot = diff(p, t);
p_ddot = diff(p_dot, t);

t = 0;
eqn1 = subs(p) - x_0;
eqn2 = subs(p_dot) - v_0;
eqn3 = subs(p_ddot) - a_0;

t = ts;
eqn4 = subs(p) - x_1;
eqn5 = subs(p_dot) - v_1;
eqn6 = subs(p_ddot) - a_1;

sol = solve([eqn1, eqn2, eqn3, eqn4, eqn5, eqn6], [x1,x2,x3,x4,x5,x6]);