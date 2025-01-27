function theta = inverse_kin(pos)
% Note: there are 8 possible solutions, 2 for every angle theta_i (one 
% "elbow up" and the other "elbow down") but we should take the solution 
% which satisfies the robot constraints (generally the solution "elbow down")

syms theta1 theta2 theta3 real;
x = pos(1);
y = pos(2);
z = pos(3);

% Robot geometric parameters
robot_geom_param

% Fictional knees coordinates in frame {R}
Pc1 = Rz(alpha1) * [R+l_A*cos(theta1), 0, -l_A*sin(theta1)].';
Pc2 = Rz(alpha2) * [R+l_A*cos(theta2), 0, -l_A*sin(theta2)].';
Pc3 = Rz(alpha3) * [R+l_A*cos(theta3), 0, -l_A*sin(theta3)].';

% Foot coordinates in frame {R}
P = [x, y, z].';
% % Ankles cohordinates in frame {R}
% Pb1 = P + Rb*[cos(alpha1), sin(alpha1), 0].';
% Pb2 = P + Rb*[cos(alpha2), sin(alpha2), 0].';
% Pb3 = P + Rb*[cos(alpha3), sin(alpha3), 0].';

% Equations || Pci-P ||^2 = l_b^2
eqn1 = (Pc1(1)-P(1))^2 + (Pc1(2)-P(2))^2 + (Pc1(3)-P(3))^2 == l_B^2;
eqn2 = (Pc2(1)-P(1))^2 + (Pc2(2)-P(2))^2 + (Pc2(3)-P(3))^2 == l_B^2;
eqn3 = (Pc3(1)-P(1))^2 + (Pc3(2)-P(2))^2 + (Pc3(3)-P(3))^2 == l_B^2;

% Solutions
th1 = solve(eqn1, theta1, 'Real', true);
th2 = solve(eqn2, theta2, 'Real', true);
th3 = solve(eqn3, theta3, 'Real', true);

% TO DO: Choose the correct solution

theta = [th1, th2, th3];