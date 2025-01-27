%% Initialization 
clear; clc;

tic
fprintf('Starting preliminary computations...\n')

addpath("./dynamics/")
addpath("./kinematics/")
addpath("./utility/")
addpath("./init_parameters/")

syms bx by bz phi zeta psi ...
     th11 th12 th13 p1x p1y p1z ...
     th21 th22 th23 p2x p2y p2z ...
     th31 th32 th33 p3x p3y p3z ...
     th41 th42 th43 p4x p4y p4z real;
syms bx_dot by_dot bz_dot phi_dot zeta_dot psi_dot ...
     th11_dot th12_dot th13_dot p1x_dot p1y_dot p1z_dot  ...
     th21_dot th22_dot th23_dot p2x_dot p2y_dot p2z_dot ...
     th31_dot th32_dot th33_dot p3x_dot p3y_dot p3z_dot ...
     th41_dot th42_dot th43_dot p4x_dot p4y_dot p4z_dot real;

q = [bx;by;bz;phi;zeta;psi;
     th11;th12;th13;p1x;p1y;p1z;
     th21;th22;th23;p2x;p2y;p2z;
     th31;th32;th33;p3x;p3y;p3z;
     th41;th42;th43;p4x;p4y;p4z];

q_dot = [bx_dot;by_dot;bz_dot;phi_dot;zeta_dot;psi_dot; ... 
         th11_dot;th12_dot;th13_dot;p1x_dot;p1y_dot;p1z_dot;
         th21_dot;th22_dot;th23_dot;p2x_dot;p2y_dot;p2z_dot;
         th31_dot;th32_dot;th33_dot;p3x_dot;p3y_dot;p3z_dot;
         th41_dot;th42_dot;th43_dot;p4x_dot;p4y_dot;p4z_dot];

g = [0; 0; -9.8]; % gravity acceleration in frame {S}

% Robot geometric parameters
robot_geom_param

% Robot dynamics parameters
deltino_dyn_param

% Upper legs inertiae w.r.t. CoM in {T_ij} frame
r_CM = (l_A/2*m_t + l_A*m_k)/(m_t + m_k);

I_thigh = [0, 0, 0;
           0, 1, 0;
           0, 0, 1]*1/12*m_t*l_A^2;

I_steiner = [0                  0                  0; 
             0 m_t*(r_CM-l_A/2)^2                  0; 
             0                  0 m_t*(r_CM-l_A/2)^2];

I_add = [0                0                0; 
         0 m_k*(l_A-r_CM)^2                0; 
         0                0 m_k*(l_A-r_CM)^2];

I_leg_body = I_thigh + I_steiner + I_add; 
 
% Trunk inertia w.r.t. CoM in {B} frame
I_trunk_body = [l2^2*h^2, 0, 0;
                0, l1^2*h^2, 0;
                0, 0, l2^2*l1^2]*1/12*m_b;

duration = toc;
fprintf('... duration = %d\n\n', duration)

%% Kinematics 
tic
fprintf('Starting computation of knematics ...\n')

% Orientation of frames {B_ij} w.r.t global frame {B}
alpha1 = 0;
alpha2 = 2/3*pi;
alpha3 = 4/3*pi; 

% Transformation matrices from {S} to {B} and {B_i}
R_SB = Rypr(phi, zeta, psi);
g_SB = [R_SB  [bx;by;bz];
        0 0 0         1 ];

R_SBi1 = R_SB*Rz(alpha1);
R_SBi2 = R_SB*Rz(alpha2);
R_SBi3 = R_SB*Rz(alpha3);

offset1 = [l1/2; l2/2; -h/2];
d_SB1 = [bx;by;bz] + R_SB*offset1;
g_SB11 = [R_SBi1  d_SB1;
          0 0 0      1];
g_SB12 = [R_SBi2  d_SB1;
          0 0 0      1];
g_SB13 = [R_SBi3  d_SB1;
          0 0 0      1];

offset2 = [l1/2; -l2/2; -h/2];
d_SB2 = [bx;by;bz] + R_SB*offset2;
g_SB21 = [R_SBi1  d_SB2;
          0 0 0      1];
g_SB22 = [R_SBi2  d_SB2;
          0 0 0      1];
g_SB23 = [R_SBi3  d_SB2;
          0 0 0      1];

offset3 = [-l1/2; -l2/2; -h/2];
d_SB3 = [bx;by;bz] + R_SB*offset3;
g_SB31 = [R_SBi1  d_SB3;
          0 0 0      1];
g_SB32 = [R_SBi2  d_SB3;
          0 0 0      1];
g_SB33 = [R_SBi3  d_SB3;
          0 0 0      1];

offset4 = [-l1/2; l2/2; -h/2];
d_SB4 = [bx;by;bz] + R_SB*offset4;
g_SB41 = [R_SBi1  d_SB4;
          0 0 0      1];
g_SB42 = [R_SBi2  d_SB4;
          0 0 0      1];
g_SB43 = [R_SBi3  d_SB4;
          0 0 0      1];

% Rotational matrices between {T_ij} frames and {S} frame
% First leg
R_S_T11 = R_SBi1*Ry(th11);
R_S_T12 = R_SBi2*Ry(th12);
R_S_T13 = R_SBi3*Ry(th13);
% Second leg
R_S_T21 = R_SBi1*Ry(th21);
R_S_T22 = R_SBi2*Ry(th22);
R_S_T23 = R_SBi3*Ry(th23);
% Third leg
R_S_T31 = R_SBi1*Ry(th31);
R_S_T32 = R_SBi2*Ry(th32);
R_S_T33 = R_SBi3*Ry(th33);
% Fourth leg
R_S_T41 = R_SBi1*Ry(th41);
R_S_T42 = R_SBi2*Ry(th42);
R_S_T43 = R_SBi3*Ry(th43);

% Upper legs inertiae w.r.t. CoM in {S} frame
% First leg
I_leg11_spatial = R_S_T11 * I_leg_body * R_S_T11'; 
I_leg12_spatial = R_S_T12 * I_leg_body * R_S_T12'; 
I_leg13_spatial = R_S_T13 * I_leg_body * R_S_T13'; 
% Second leg
I_leg21_spatial = R_S_T21 * I_leg_body * R_S_T21'; 
I_leg22_spatial = R_S_T22 * I_leg_body * R_S_T22'; 
I_leg23_spatial = R_S_T23 * I_leg_body * R_S_T23'; 
% Third leg
I_leg31_spatial = R_S_T31 * I_leg_body * R_S_T31'; 
I_leg32_spatial = R_S_T32 * I_leg_body * R_S_T32'; 
I_leg33_spatial = R_S_T33 * I_leg_body * R_S_T33'; 
% Fourth leg
I_leg41_spatial = R_S_T41 * I_leg_body * R_S_T41'; 
I_leg42_spatial = R_S_T42 * I_leg_body * R_S_T42'; 
I_leg43_spatial = R_S_T43 * I_leg_body * R_S_T43'; 

% Trunk inertia w.r.t. CoM in {S} frame 
I_trunk_spatial = R_SB * I_trunk_body * R_SB';

duration = toc;
fprintf('... duration = %d\n\n', duration)

%% Differential Kinematics

tic
fprintf('Starting computation of differential kinematics ...\n\n')

% Jacobians
% Central body (Trunk)
J_vb = [eye(3), zeros(3, 27)];
J_wb = [zeros(3), [0;0;1], Rz(phi)*[0;1;0], Rz(phi)*Ry(zeta)*[1;0;0], ...
        zeros(3, 24)];
J_b = [J_vb; J_wb];

% First delta
% -- First thigh
d_B11T11 = [R_A + r_CM*cos(th11); 0; -r_CM*sin(th11)];
g_B11T11 = [Ry(th11), d_B11T11;
             0, 0, 0,        1];
g_ST11 = g_SB11*g_B11T11;
d_ST11 = g_ST11(1:3,4);
J_v11 = simplify(jacobian(d_ST11, q), 'steps', 10);
J_w11 = [zeros(3), [0;0;1], Rz(phi)*[0;1;0], Rz(phi)*Ry(zeta)*[1;0;0], ...
        R_SBi1*[0;1;0], zeros(3,23)];
J11 = [J_v11; J_w11];
% -- Second thigh
d_B12T12 = [R_A + r_CM*cos(th12); 0; -r_CM*sin(th12)];
g_B12T12 = [Ry(th12), d_B12T12;
             0, 0, 0,        1];
g_ST12 = g_SB12*g_B12T12;
d_ST12 = g_ST12(1:3,4);
J_v12 = simplify(jacobian(d_ST12, q), 'steps', 10);
J_w12 = [zeros(3), [0;0;1], Rz(phi)*[0;1;0], Rz(phi)*Ry(zeta)*[1;0;0], ...
        zeros(3,1), R_SBi2*[0;1;0], zeros(3,22)];
J12 = [J_v12; J_w12];
% -- Third thigh
d_B13T13 = [R_A + r_CM*cos(th13); 0; -r_CM*sin(th13)];
g_B13T13 = [Ry(th13), d_B13T13;
             0, 0, 0,        1];
g_ST13 = g_SB13*g_B13T13;
d_ST13 = g_ST13(1:3,4);
J_v13 = simplify(jacobian(d_ST13, q), 'steps', 10);
J_w13 = [zeros(3), [0;0;1], Rz(phi)*[0;1;0], Rz(phi)*Ry(zeta)*[1;0;0], ...
        zeros(3,2), R_SBi3*[0;1;0], zeros(3,21)];
J13 = [J_v13; J_w13];
% -- Coupler
J_v14 = [zeros(3,9), eye(3), zeros(3,18)];
J_w14 = [zeros(3), [0;0;1], Rz(phi)*[0;1;0], Rz(phi)*Ry(zeta)*[1;0;0], ...
        zeros(3, 24)];
J14 = [J_v14; J_w14];

% Second delta
% -- First thigh
d_B21T21 = [R_A + r_CM*cos(th21); 0; -r_CM*sin(th21)];
g_B21T21 = [Ry(th21), d_B21T21;
             0, 0, 0,        1];
g_ST21 = g_SB21*g_B21T21;
d_ST21 = g_ST21(1:3,4);
J_v21 = simplify(jacobian(d_ST21, q), 'steps', 10);
J_w21 = [zeros(3), [0;0;1], Rz(phi)*[0;1;0], Rz(phi)*Ry(zeta)*[1;0;0], ...
         zeros(3,6), R_SBi1*[0;1;0], zeros(3,17)];
J21 = [J_v21; J_w21];
% -- Second thigh
d_B22T22 = [R_A + r_CM*cos(th22); 0; -r_CM*sin(th22)];
g_B22T22 = [Ry(th22), d_B22T22;
             0, 0, 0,        1];
g_ST22 = g_SB22*g_B22T22;
d_ST22 = g_ST22(1:3,4);
J_v22 = simplify(jacobian(d_ST22, q), 'steps', 10);
J_w22 = [zeros(3), [0;0;1], Rz(phi)*[0;1;0], Rz(phi)*Ry(zeta)*[1;0;0], ...
         zeros(3,7), R_SBi2*[0;1;0], zeros(3,16)];
J22 = [J_v22; J_w22];
% -- Third thigh
d_B23T23 = [R_A + r_CM*cos(th23); 0; -r_CM*sin(th23)];
g_B23T23 = [Ry(th23), d_B23T23;
             0, 0, 0,        1];
g_ST23 = g_SB23*g_B23T23;
d_ST23 = g_ST23(1:3,4);
J_v23 = simplify(jacobian(d_ST23, q), 'steps', 10);
J_w23 = [zeros(3), [0;0;1], Rz(phi)*[0;1;0], Rz(phi)*Ry(zeta)*[1;0;0], ...
        zeros(3,8), R_SBi3*[0;1;0], zeros(3,15)];
J23 = [J_v23; J_w23];
% -- Coupler
J_v24 = [zeros(3,15), eye(3), zeros(3,12)];
J_w24 = [zeros(3), [0;0;1], Rz(phi)*[0;1;0], Rz(phi)*Ry(zeta)*[1;0;0], ...
        zeros(3, 24)];
J24 = [J_v24; J_w24];

% Third delta
% -- First thigh
d_B31T31 = [R_A + r_CM*cos(th31); 0; -r_CM*sin(th31)];
g_B31T31 = [Ry(th31), d_B31T31;
             0, 0, 0,        1];
g_ST31 = g_SB31*g_B31T31;
d_ST31 = g_ST31(1:3,4);
J_v31 = simplify(jacobian(d_ST31, q), 'steps', 10);
J_w31 = [zeros(3), [0;0;1], Rz(phi)*[0;1;0], Rz(phi)*Ry(zeta)*[1;0;0], ...
         zeros(3,12), R_SBi1*[0;1;0], zeros(3,11)];
J31 = [J_v31; J_w31];
% -- Second thigh
d_B32T32 = [R_A + r_CM*cos(th32); 0; -r_CM*sin(th32)];
g_B32T32 = [Ry(th32), d_B32T32;
             0, 0, 0,        1];
g_ST32 = g_SB32*g_B32T32;
d_ST32 = g_ST32(1:3,4);
J_v32 = simplify(jacobian(d_ST32, q), 'steps', 10);
J_w32 = [zeros(3), [0;0;1], Rz(phi)*[0;1;0], Rz(phi)*Ry(zeta)*[1;0;0], ...
         zeros(3,13), R_SBi2*[0;1;0], zeros(3,10)];
J32 = [J_v32; J_w32];
% -- Third thigh
d_B33T33 = [R_A + r_CM*cos(th33); 0; -r_CM*sin(th33)];
g_B33T33 = [Ry(th33), d_B33T33;
             0, 0, 0,        1];
g_ST33 = g_SB33*g_B33T33;
d_ST33 = g_ST33(1:3,4);
J_v33 = simplify(jacobian(d_ST33, q), 'steps', 10);
J_w33 = [zeros(3), [0;0;1], Rz(phi)*[0;1;0], Rz(phi)*Ry(zeta)*[1;0;0], ...
        zeros(3,14), R_SBi3*[0;1;0], zeros(3,9)];
J33 = [J_v33; J_w33];
% -- Coupler
J_v34 = [zeros(3,21), eye(3), zeros(3,6)];
J_w34 = [zeros(3), [0;0;1], Rz(phi)*[0;1;0], Rz(phi)*Ry(zeta)*[1;0;0], ...
        zeros(3, 24)];
J34 = [J_v34; J_w34];

% Fourth delta 
% -- First thigh
d_B41T41 = [R_A + r_CM*cos(th41); 0; -r_CM*sin(th41)];
g_B41T41 = [Ry(th41), d_B41T41;
             0, 0, 0,        1];
g_ST41 = g_SB41*g_B41T41;
d_ST41 = g_ST41(1:3,4);
J_v41 = simplify(jacobian(d_ST41, q), 'steps', 10);
J_w41 = [zeros(3), [0;0;1], Rz(phi)*[0;1;0], Rz(phi)*Ry(zeta)*[1;0;0], ...
         zeros(3,18), R_SBi1*[0;1;0], zeros(3,5)];
J41 = [J_v41; J_w41];
% -- Second thigh
d_B42T42 = [R_A + r_CM*cos(th42); 0; -r_CM*sin(th42)];
g_B42T42 = [Ry(th42), d_B42T42;
             0, 0, 0,        1];
g_ST42 = g_SB42*g_B42T42;
d_ST42 = g_ST42(1:3,4);
J_v42 = simplify(jacobian(d_ST42, q), 'steps', 10);
J_w42 = [zeros(3), [0;0;1], Rz(phi)*[0;1;0], Rz(phi)*Ry(zeta)*[1;0;0], ...
         zeros(3,19), R_SBi2*[0;1;0], zeros(3,4)];
J42 = [J_v42; J_w42];
% -- Third thigh
d_B43T43 = [R_A + r_CM*cos(th43); 0; -r_CM*sin(th43)];
g_B43T43 = [Ry(th43), d_B43T43;
             0, 0, 0,        1];
g_ST43 = g_SB43*g_B43T43;
d_ST43 = g_ST43(1:3,4);
J_v43 = simplify(jacobian(d_ST43, q), 'steps', 10);
J_w43 = [zeros(3), [0;0;1], Rz(phi)*[0;1;0], Rz(phi)*Ry(zeta)*[1;0;0], ...
        zeros(3,20), R_SBi3*[0;1;0], zeros(3,3)];
J43 = [J_v43; J_w43];
% -- Coupler
J_v44 = [zeros(3,27), eye(3)];
J_w44 = [zeros(3), [0;0;1], Rz(phi)*[0;1;0], Rz(phi)*Ry(zeta)*[1;0;0], ...
        zeros(3, 24)];
J44 = [J_v44; J_w44];

duration = toc;
fprintf('... duration = %d\n\n', duration)

%% Dynamics
fprintf('Starting computation of dynamical matrices...\n')
tic
% Compute dynamic matrices
jacobians = [J_b, J11, J12, J13, J14, ...
             J21, J22, J23, J24, ...
             J31, J32, J33, J34, ...
             J41, J42, J43, J44];
jacobians = simplify(jacobians, 'steps', 100);
masses = [m_b, m_l, m_l, m_l, m_c, ...
          m_l, m_l, m_l, m_c,...
          m_l, m_l, m_l, m_c,...
          m_l, m_l, m_l, m_c];
inertiae = [I_trunk_spatial, ...
            I_leg11_spatial, I_leg12_spatial, I_leg13_spatial, zeros(3),...
            I_leg21_spatial, I_leg22_spatial, I_leg23_spatial, zeros(3),...
            I_leg31_spatial, I_leg32_spatial, I_leg33_spatial, zeros(3),...
            I_leg41_spatial, I_leg42_spatial, I_leg43_spatial, zeros(3)];
        
tic
fprintf('Starting computation of M matrix...\n')
M = compute_M_parallel(masses, inertiae, jacobians);
duration = toc;
fprintf('... duration = %d\n\n', duration)

tic
fprintf('Starting computation of C matrix...\n')
C = compute_C_parallel(M, q, q_dot);
duration = toc;
fprintf('... duration = %d\n\n', duration)

tic
fprintf('Starting computation of G vector...\n')
G = compute_G_parallel(masses, jacobians, g);
duration = toc;
fprintf('... duration = %d\n\n', duration)

%% Constraints
tic
fprintf('Starting computation of A matrix...\n')

% Compute constraints matrix
% First delta
s11 = [p1x;p1y;p1z;1] - g_SB11*[R+l_A*cos(th11); 0; -l_A*sin(th11); 1];
s12 = [p1x;p1y;p1z;1] - g_SB12*[R+l_A*cos(th12); 0; -l_A*sin(th12); 1];
s13 = [p1x;p1y;p1z;1] - g_SB13*[R+l_A*cos(th13); 0; -l_A*sin(th13); 1];

v11 = [R+l_A*cos(th11); 0; -l_A*sin(th11)];
v12 = [R+l_A*cos(th12); 0; -l_A*sin(th12)];
v13 = [R+l_A*cos(th13); 0; -l_A*sin(th13)];

c11 = [diff(Rz(phi))*Ry(zeta)*Rx(psi)*Rz(alpha1)*v11 + ...
       diff(Rz(phi))*Ry(zeta)*Rx(psi)*offset1; 0];
c12 = [diff(Rz(phi))*Ry(zeta)*Rx(psi)*Rz(alpha2)*v12 + ...
       diff(Rz(phi))*Ry(zeta)*Rx(psi)*offset1; 0];
c13 = [diff(Rz(phi))*Ry(zeta)*Rx(psi)*Rz(alpha3)*v13 + ...
       diff(Rz(phi))*Ry(zeta)*Rx(psi)*offset1; 0];

d11 = [Rz(phi)*diff(Ry(zeta))*Rx(psi)*Rz(alpha1)*v11 + ...
       Rz(phi)*diff(Ry(zeta))*Rx(psi)*offset1; 0];
d12 = [Rz(phi)*diff(Ry(zeta))*Rx(psi)*Rz(alpha2)*v12 + ...
       Rz(phi)*diff(Ry(zeta))*Rx(psi)*offset1; 0];
d13 = [Rz(phi)*diff(Ry(zeta))*Rx(psi)*Rz(alpha3)*v13 + ...
       Rz(phi)*diff(Ry(zeta))*Rx(psi)*offset1; 0];

e11 = [Rz(phi)*Ry(zeta)*diff(Rx(psi))*Rz(alpha1)*v11 + ...
       Rz(phi)*Ry(zeta)*diff(Rx(psi))*offset1; 0];
e12 = [Rz(phi)*Ry(zeta)*diff(Rx(psi))*Rz(alpha2)*v12 + ...
       Rz(phi)*Ry(zeta)*diff(Rx(psi))*offset1; 0];
e13 = [Rz(phi)*Ry(zeta)*diff(Rx(psi))*Rz(alpha3)*v13 + ...
       Rz(phi)*Ry(zeta)*diff(Rx(psi))*offset1; 0];

f11 = g_SB11 * [-l_A*sin(th11); 0; -l_A*cos(th11); 0];
f12 = g_SB12 * [-l_A*sin(th12); 0; -l_A*cos(th12); 0];
f13 = g_SB13 * [-l_A*sin(th13); 0; -l_A*cos(th13); 0];

A1 = [s11(1) s11(2) s11(3) s11'*c11 s11'*d11 s11'*e11 s11'*f11 0 0 -s11(1) -s11(2) -s11(3), zeros(1,18);
      s12(1) s12(2) s12(3) s12'*c12 s12'*d12 s12'*e12 0 s12'*f12 0 -s12(1) -s12(2) -s12(3), zeros(1,18);
      s13(1) s13(2) s13(3) s13'*c13 s13'*d13 s13'*e13 0 0 s13'*f13 -s13(1) -s13(2) -s13(3), zeros(1,18)];

% Second delta
s21 = [p2x;p2y;p2z;1] - g_SB21*[R+l_A*cos(th21); 0; -l_A*sin(th21); 1];
s22 = [p2x;p2y;p2z;1] - g_SB22*[R+l_A*cos(th22); 0; -l_A*sin(th22); 1];
s23 = [p2x;p2y;p2z;1] - g_SB23*[R+l_A*cos(th23); 0; -l_A*sin(th23); 1];

v21 = [R+l_A*cos(th21); 0; -l_A*sin(th21)];
v22 = [R+l_A*cos(th22); 0; -l_A*sin(th22)];
v23 = [R+l_A*cos(th23); 0; -l_A*sin(th23)];

c21 = [diff(Rz(phi))*Ry(zeta)*Rx(psi)*Rz(alpha1)*v21 + ...
       diff(Rz(phi))*Ry(zeta)*Rx(psi)*offset2; 0];
c22 = [diff(Rz(phi))*Ry(zeta)*Rx(psi)*Rz(alpha2)*v22 + ...
       diff(Rz(phi))*Ry(zeta)*Rx(psi)*offset2; 0];
c23 = [diff(Rz(phi))*Ry(zeta)*Rx(psi)*Rz(alpha3)*v23 + ...
       diff(Rz(phi))*Ry(zeta)*Rx(psi)*offset2; 0];

d21 = [Rz(phi)*diff(Ry(zeta))*Rx(psi)*Rz(alpha1)*v21 + ...
       Rz(phi)*diff(Ry(zeta))*Rx(psi)*offset2; 0];
d22 = [Rz(phi)*diff(Ry(zeta))*Rx(psi)*Rz(alpha2)*v22 + ...
       Rz(phi)*diff(Ry(zeta))*Rx(psi)*offset2; 0];
d23 = [Rz(phi)*diff(Ry(zeta))*Rx(psi)*Rz(alpha3)*v23 + ...
       Rz(phi)*diff(Ry(zeta))*Rx(psi)*offset2; 0];

e21 = [Rz(phi)*Ry(zeta)*diff(Rx(psi))*Rz(alpha1)*v21 + ...
       Rz(phi)*Ry(zeta)*diff(Rx(psi))*offset2; 0];
e22 = [Rz(phi)*Ry(zeta)*diff(Rx(psi))*Rz(alpha2)*v22 + ...
       Rz(phi)*Ry(zeta)*diff(Rx(psi))*offset2; 0];
e23 = [Rz(phi)*Ry(zeta)*diff(Rx(psi))*Rz(alpha3)*v23 + ...
       Rz(phi)*Ry(zeta)*diff(Rx(psi))*offset2; 0];

f21 = g_SB21 * [-l_A*sin(th21); 0; -l_A*cos(th21); 0];
f22 = g_SB22 * [-l_A*sin(th22); 0; -l_A*cos(th22); 0];
f23 = g_SB23 * [-l_A*sin(th23); 0; -l_A*cos(th23); 0];

A2 = [s21(1) s21(2) s21(3) s21'*c21 s21'*d21 s21'*e21 zeros(1,6) s21'*f21 0 0 -s21(1) -s21(2) -s21(3), zeros(1, 12);
      s22(1) s22(2) s22(3) s22'*c22 s22'*d22 s22'*e22 zeros(1,6) 0 s22'*f22 0 -s22(1) -s22(2) -s22(3), zeros(1, 12);
      s23(1) s23(2) s23(3) s23'*c23 s23'*d23 s23'*e23 zeros(1,6) 0 0 s23'*f23 -s23(1) -s23(2) -s23(3), zeros(1, 12)];

% Third delta
s31 = [p3x;p3y;p3z;1] - g_SB31*[R+l_A*cos(th31); 0; -l_A*sin(th31); 1];
s32 = [p3x;p3y;p3z;1] - g_SB32*[R+l_A*cos(th32); 0; -l_A*sin(th32); 1];
s33 = [p3x;p3y;p3z;1] - g_SB33*[R+l_A*cos(th33); 0; -l_A*sin(th33); 1];

v31 = [R+l_A*cos(th31); 0; -l_A*sin(th31)];
v32 = [R+l_A*cos(th32); 0; -l_A*sin(th32)];
v33 = [R+l_A*cos(th33); 0; -l_A*sin(th33)];

c31 = [diff(Rz(phi))*Ry(zeta)*Rx(psi)*Rz(alpha1)*v31 + ...
       diff(Rz(phi))*Ry(zeta)*Rx(psi)*offset3; 0];
c32 = [diff(Rz(phi))*Ry(zeta)*Rx(psi)*Rz(alpha2)*v32 + ...
       diff(Rz(phi))*Ry(zeta)*Rx(psi)*offset3; 0];
c33 = [diff(Rz(phi))*Ry(zeta)*Rx(psi)*Rz(alpha3)*v33 + ...
       diff(Rz(phi))*Ry(zeta)*Rx(psi)*offset3; 0];

d31 = [Rz(phi)*diff(Ry(zeta))*Rx(psi)*Rz(alpha1)*v31 + ...
       Rz(phi)*diff(Ry(zeta))*Rx(psi)*offset3; 0];
d32 = [Rz(phi)*diff(Ry(zeta))*Rx(psi)*Rz(alpha2)*v32 + ...
       Rz(phi)*diff(Ry(zeta))*Rx(psi)*offset3; 0];
d33 = [Rz(phi)*diff(Ry(zeta))*Rx(psi)*Rz(alpha3)*v33 + ...
       Rz(phi)*diff(Ry(zeta))*Rx(psi)*offset3; 0];

e31 = [Rz(phi)*Ry(zeta)*diff(Rx(psi))*Rz(alpha1)*v31 + ...
       Rz(phi)*Ry(zeta)*diff(Rx(psi))*offset3; 0];
e32 = [Rz(phi)*Ry(zeta)*diff(Rx(psi))*Rz(alpha2)*v32 + ...
       Rz(phi)*Ry(zeta)*diff(Rx(psi))*offset3; 0];
e33 = [Rz(phi)*Ry(zeta)*diff(Rx(psi))*Rz(alpha3)*v33 + ...
       Rz(phi)*Ry(zeta)*diff(Rx(psi))*offset3; 0];

f31 = g_SB31 * [-l_A*sin(th31); 0; -l_A*cos(th31); 0];
f32 = g_SB32 * [-l_A*sin(th32); 0; -l_A*cos(th32); 0];
f33 = g_SB33 * [-l_A*sin(th33); 0; -l_A*cos(th33); 0];

A3 = [s31(1) s31(2) s31(3) s31'*c31 s31'*d31 s31'*e31 zeros(1,12) s31'*f31 0 0 -s31(1) -s31(2) -s31(3), zeros(1,6);
      s32(1) s32(2) s32(3) s32'*c32 s32'*d32 s32'*e32 zeros(1,12) 0 s32'*f32 0 -s32(1) -s32(2) -s32(3), zeros(1,6);
      s33(1) s33(2) s33(3) s33'*c33 s33'*d33 s33'*e33 zeros(1,12) 0 0 s33'*f33 -s33(1) -s33(2) -s33(3), zeros(1,6)];

% Fourth delta
s41 = [p4x;p4y;p4z;1] - g_SB41*[R+l_A*cos(th41); 0; -l_A*sin(th41); 1];
s42 = [p4x;p4y;p4z;1] - g_SB42*[R+l_A*cos(th42); 0; -l_A*sin(th42); 1];
s43 = [p4x;p4y;p4z;1] - g_SB43*[R+l_A*cos(th43); 0; -l_A*sin(th43); 1];

v41 = [R+l_A*cos(th41); 0; -l_A*sin(th41)];
v42 = [R+l_A*cos(th42); 0; -l_A*sin(th42)];
v43 = [R+l_A*cos(th43); 0; -l_A*sin(th43)];

c41 = [diff(Rz(phi))*Ry(zeta)*Rx(psi)*Rz(alpha1)*v41 + ...
       diff(Rz(phi))*Ry(zeta)*Rx(psi)*offset4; 0];
c42 = [diff(Rz(phi))*Ry(zeta)*Rx(psi)*Rz(alpha2)*v42 + ...
       diff(Rz(phi))*Ry(zeta)*Rx(psi)*offset4; 0];
c43 = [diff(Rz(phi))*Ry(zeta)*Rx(psi)*Rz(alpha3)*v43 + ...
       diff(Rz(phi))*Ry(zeta)*Rx(psi)*offset4; 0];

d41 = [Rz(phi)*diff(Ry(zeta))*Rx(psi)*Rz(alpha1)*v41 + ...
       Rz(phi)*diff(Ry(zeta))*Rx(psi)*offset4; 0];
d42 = [Rz(phi)*diff(Ry(zeta))*Rx(psi)*Rz(alpha2)*v42 + ...
       Rz(phi)*diff(Ry(zeta))*Rx(psi)*offset4; 0];
d43 = [Rz(phi)*diff(Ry(zeta))*Rx(psi)*Rz(alpha3)*v43 + ...
       Rz(phi)*diff(Ry(zeta))*Rx(psi)*offset4; 0];

e41 = [Rz(phi)*Ry(zeta)*diff(Rx(psi))*Rz(alpha1)*v41 + ...
       Rz(phi)*Ry(zeta)*diff(Rx(psi))*offset4; 0];
e42 = [Rz(phi)*Ry(zeta)*diff(Rx(psi))*Rz(alpha2)*v42 + ...
       Rz(phi)*Ry(zeta)*diff(Rx(psi))*offset4; 0];
e43 = [Rz(phi)*Ry(zeta)*diff(Rx(psi))*Rz(alpha3)*v43 + ...
       Rz(phi)*Ry(zeta)*diff(Rx(psi))*offset4; 0];

f41 = g_SB41 * [-l_A*sin(th41); 0; -l_A*cos(th41); 0];
f42 = g_SB42 * [-l_A*sin(th42); 0; -l_A*cos(th42); 0];
f43 = g_SB43 * [-l_A*sin(th43); 0; -l_A*cos(th43); 0];

A4 = [s41(1) s41(2) s41(3) s41'*c41 s41'*d41 s41'*e41 zeros(1,18) s41'*f41 0 0 -s41(1) -s41(2) -s41(3);
      s42(1) s42(2) s42(3) s42'*c42 s42'*d42 s42'*e42 zeros(1,18) 0 s42'*f42 0 -s42(1) -s42(2) -s42(3);
      s43(1) s43(2) s43(3) s43'*c43 s43'*d43 s43'*e43 zeros(1,18) 0 0 s43'*f43 -s43(1) -s43(2) -s43(3)];

% Final constraint matrix
A = [A1; A2; A3; A4];
duration = toc;
fprintf('... duration = %d\n\n', duration)

%% Compute a basis for nullspace of A ...

tic
fprintf('Starting computation of S matrix...\n')

S = [-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 
     0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 
     0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 
     0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 
     0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0; 
     0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0; 
     s11(1)/(s11'*f11) s11(2)/(s11'*f11) s11(3)/(s11'*f11) s11'*c11/(s11'*f11) s11'*d11/(s11'*f11) s11'*e11/(s11'*f11) s11(1)/(s11'*f11) s11(2)/(s11'*f11) s11(3)/(s11'*f11) 0 0 0 0 0 0 0 0 0;
     s12(1)/(s12'*f12) s12(2)/(s12'*f12) s12(3)/(s12'*f12) s12'*c12/(s12'*f12) s12'*d12/(s12'*f12) s12'*e12/(s12'*f12) s12(1)/(s12'*f12) s12(2)/(s12'*f12) s12(3)/(s12'*f12) 0 0 0 0 0 0 0 0 0;
     s13(1)/(s13'*f13) s13(2)/(s13'*f13) s13(3)/(s13'*f13) s13'*c13/(s13'*f13) s13'*d13/(s13'*f13) s13'*e13/(s13'*f13) s13(1)/(s13'*f13) s13(2)/(s13'*f13) s13(3)/(s13'*f13) 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;
     s21(1)/(s21'*f21) s21(2)/(s21'*f21) s21(3)/(s21'*f21) s21'*c21/(s21'*f21) s21'*d21/(s21'*f21) s21'*e21/(s21'*f21) 0 0 0 s21(1)/(s21'*f21) s21(2)/(s21'*f21) s21(3)/(s21'*f21) 0 0 0 0 0 0;
     s22(1)/(s22'*f22) s22(2)/(s22'*f22) s22(3)/(s22'*f22) s22'*c22/(s22'*f22) s22'*d22/(s22'*f22) s22'*e22/(s22'*f22) 0 0 0 s22(1)/(s22'*f22) s22(2)/(s22'*f22) s22(3)/(s22'*f22) 0 0 0 0 0 0;
     s23(1)/(s23'*f23) s23(2)/(s23'*f23) s23(3)/(s23'*f23) s23'*c23/(s23'*f23) s23'*d23/(s23'*f23) s23'*e23/(s23'*f23) 0 0 0 s23(1)/(s23'*f23) s23(2)/(s23'*f23) s23(3)/(s23'*f23) 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0;
     s31(1)/(s31'*f31) s31(2)/(s31'*f31) s31(3)/(s31'*f31) s31'*c31/(s31'*f31) s31'*d31/(s31'*f31) s31'*e31/(s31'*f31) 0 0 0 0 0 0 s31(1)/(s31'*f31) s31(2)/(s31'*f31) s31(3)/(s31'*f31) 0 0 0;
     s32(1)/(s32'*f32) s32(2)/(s32'*f32) s32(3)/(s32'*f32) s32'*c32/(s32'*f32) s32'*d32/(s32'*f32) s32'*e32/(s32'*f32) 0 0 0 0 0 0 s32(1)/(s32'*f32) s32(2)/(s32'*f32) s32(3)/(s32'*f32) 0 0 0;
     s33(1)/(s33'*f33) s33(2)/(s33'*f33) s33(3)/(s33'*f33) s33'*c33/(s33'*f33) s33'*d33/(s33'*f33) s33'*e33/(s33'*f33) 0 0 0 0 0 0 s33(1)/(s33'*f33) s33(2)/(s33'*f33) s33(3)/(s33'*f33) 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0;
     s41(1)/(s41'*f41) s41(2)/(s41'*f41) s41(3)/(s41'*f41) s41'*c41/(s41'*f41) s41'*d41/(s41'*f41) s41'*e41/(s41'*f41) 0 0 0 0 0 0 0 0 0 s41(1)/(s41'*f41) s41(2)/(s41'*f41) s41(3)/(s41'*f41);
     s42(1)/(s42'*f42) s42(2)/(s42'*f42) s42(3)/(s42'*f42) s42'*c42/(s42'*f42) s42'*d42/(s42'*f42) s42'*e42/(s42'*f42) 0 0 0 0 0 0 0 0 0 s42(1)/(s42'*f42) s42(2)/(s42'*f42) s42(3)/(s42'*f42);
     s43(1)/(s43'*f43) s43(2)/(s43'*f43) s43(3)/(s43'*f43) s43'*c43/(s43'*f43) s43'*d43/(s43'*f43) s43'*e43/(s43'*f43) 0 0 0 0 0 0 0 0 0 s43(1)/(s43'*f43) s43(2)/(s43'*f43) s43(3)/(s43'*f43);
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];

duration = toc;
fprintf('... duration = %d\n\n', duration)

%% ... and its derivative

tic
fprintf('Starting computation of S_dot matrix...\n')

syms t real;
syms bx_t(t) by_t(t) bz_t(t) phi_t(t) zeta_t(t) psi_t(t) ...
     th11_t(t) th12_t(t) th13_t(t) p1x_t(t) p1y_t(t) p1z_t(t) ...
     th21_t(t) th22_t(t) th23_t(t) p2x_t(t) p2y_t(t) p2z_t(t) ... 
     th31_t(t) th32_t(t) th33_t(t) p3x_t(t) p3y_t(t) p3z_t(t) ...
     th41_t(t) th42_t(t) th43_t(t) p4x_t(t) p4y_t(t) p4z_t(t);

S_tmp = subs(S, [bx, by, bz, phi, zeta, psi, th11, th12, th13, p1x, p1y, p1z ...
                th21, th22, th23, p2x, p2y, p2z, ...
                th31, th32, th33, p3x, p3y, p3z ...
                th41, th42, th43, p4x, p4y, p4z], ...
            [bx_t(t) by_t(t) bz_t(t) phi_t(t) zeta_t(t) psi_t(t) ...
             th11_t(t) th12_t(t) th13_t(t) p1x_t(t) p1y_t(t) p1z_t(t) ...
             th21_t(t) th22_t(t) th23_t(t) p2x_t(t) p2y_t(t) p2z_t(t) ... 
             th31_t(t) th32_t(t) th33_t(t) p3x_t(t) p3y_t(t) p3z_t(t) ...
             th41_t(t) th42_t(t) th43_t(t) p4x_t(t) p4y_t(t) p4z_t(t)]);

S_dot_tmp = diff(S_tmp, t);

S_dot = subs(S_dot_tmp, ...
    [bx_t(t) by_t(t) bz_t(t) phi_t(t) zeta_t(t) psi_t(t) ...
     th11_t(t) th12_t(t) th13_t(t) p1x_t(t) p1y_t(t) p1z_t(t) ...
     th21_t(t) th22_t(t) th23_t(t) p2x_t(t) p2y_t(t) p2z_t(t) ... 
     th31_t(t) th32_t(t) th33_t(t) p3x_t(t) p3y_t(t) p3z_t(t) ...
     th41_t(t) th42_t(t) th43_t(t) p4x_t(t) p4y_t(t) p4z_t(t), ...
     diff(bx_t,t), diff(by_t,t), diff(bz_t,t), diff(phi_t,t), diff(zeta_t,t), diff(psi_t,t),...
     diff(th11_t,t), diff(th12_t,t), diff(th13_t,t), diff(p1x_t,t), diff(p1y_t,t), diff(p1z_t,t), ... 
     diff(th21_t,t), diff(th22_t,t), diff(th23_t,t), diff(p2x_t,t), diff(p2y_t,t), diff(p2z_t,t), ...
     diff(th31_t,t), diff(th32_t,t), diff(th33_t,t), diff(p3x_t,t), diff(p3y_t,t), diff(p3z_t,t), ...
     diff(th41_t,t), diff(th42_t,t), diff(th43_t,t), diff(p4x_t,t), diff(p4y_t,t), diff(p4z_t,t)], ...
    [bx, by, bz, phi, zeta, psi, ...
     th11, th12, th13, p1x, p1y, p1z ...
     th21, th22, th23, p2x, p2y, p2z, ...
     th31, th32, th33, p3x, p3y, p3z ...
     th41, th42, th43, p4x, p4y, p4z, ...
     bx_dot, by_dot, bz_dot, phi_dot, zeta_dot, psi_dot, ...
     th11_dot, th12_dot, th13_dot, p1x_dot, p1y_dot, p1z_dot, ...
     th21_dot, th22_dot, th23_dot, p2x_dot, p2y_dot, p2z_dot, ...
     th31_dot, th32_dot, th33_dot, p3x_dot, p3y_dot, p3z_dot, ...
     th41_dot, th42_dot, th43_dot, p4x_dot, p4y_dot, p4z_dot]);
 
duration = toc;
fprintf('... duration = %d\n\n', duration)

%% Simplifications
fprintf('Starting simplification of M, C, G, A, S and S_dot...\n')
temp = {};
parfor idx = 1:6
    tic
    switch idx
        case 1
            temp{idx} = simplify(M, 'steps', 100);
            duration = toc;
            fprintf('... duration for M = %d\n\n', duration)
        case 2
            temp{idx} = simplify(C, 'steps', 100);
            duration = toc;
            fprintf('... duration for C = %d\n\n', duration)
        case 3
            temp{idx} = simplify(G, 'steps', 100);
            duration = toc;
            fprintf('... duration for G= %d\n\n', duration)
        case 4
            temp{idx} = simplify(A, 'steps', 100);
            duration = toc;
            fprintf('... duration for A = %d\n\n', duration)
        case 5
            temp{idx} = simplify(S, 'steps', 100);
            duration = toc;
            fprintf('... duration for S = %d\n\n', duration)
        case 6
            temp{idx} = simplify(S_dot, 'steps', 5);
            duration = toc;
            fprintf('... duration for S_dot = %d\n\n', duration)
    end
end

M = temp{1};
C = temp{2};
G = temp{3};
A = temp{4};
S = temp{5};
S_dot = temp{6};
