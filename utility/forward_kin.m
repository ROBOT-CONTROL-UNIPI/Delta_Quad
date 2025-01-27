function pos = forward_kin(theta)
% Note: there are 2 possible solutions, one with the end effector above the
% base, and the other with the end effector below the base. We have to
% choose the l_Atter one.
addpath("kinematics");

syms x y z real;
theta1 = theta(1);
theta2 = theta(2);
theta3 = theta(3);

% Robot geometric parameters
deltino_kin_param

% Knees cohordinates in frame {R}
Pc1 = Rz(alpha1) * [R+l_A*cos(theta1), 0, -l_A*sin(theta1)].';
Pc2 = Rz(alpha2) * [R+l_A*cos(theta2), 0, -l_A*sin(theta2)].';
Pc3 = Rz(alpha3) * [R+l_A*cos(theta3), 0, -l_A*sin(theta3)].';

% Equations of tril_Ateration to be solved
eqn1 = (x - Pc1(1))^2 + (y - Pc1(2))^2 + (z - Pc1(3))^2 == l_B^2;
eqn2 = (x - Pc2(1))^2 + (y - Pc2(2))^2 + (z - Pc2(3))^2 == l_B^2;
eqn3 = (x - Pc3(1))^2 + (y - Pc3(2))^2 + (z - Pc3(3))^2 == l_B^2;
sols = solve([eqn1, eqn2, eqn3], [x, y, z]);

% Choose the correct solution
if(length(sols.z) > 1) 
    if(sols.z(1) < sols.z(2))
        px = double(sols.x(1));
        py = double(sols.y(1));
        pz = double(sols.z(1));
    else 
        px = double(sols.x(2));
        py = double(sols.y(2));
        pz = double(sols.z(2));
    end
else 
    px = double(sols.x);
    py = double(sols.y);
    pz = double(sols.z);

end
pos = [px;py;pz];