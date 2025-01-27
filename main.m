%% Path Initialization
clc; clear;
addpath("./kinematics/");
addpath("./init_parameters/");
addpath("./utility/");

%% Robot Geometric and Dynamics Parameters
deltino_kin_param
deltino_dyn_param

%% Jacobian for Contact Forces
J_1 = [zeros(3, 9), eye(3), zeros(3, 18)];
J_2 = [zeros(3, 15), eye(3), zeros(3, 12)];
J_3 = [zeros(3, 21), eye(3), zeros(3, 6)];
J_4 = [zeros(3, 27), eye(3)]; 
J_T = [J_1', J_2', J_3', J_4'];

%% Initialization of Starting Position 
% Base
base_0 = [0; 0; 0.5];
ypr_0 = [0; 0; 0];
R_SB = Rypr(ypr_0(1), ypr_0(2), ypr_0(3));
offset1 = [+l1/2; +l2/2; -h/2];
offset2 = [+l1/2; -l2/2; -h/2];
offset3 = [-l1/2; -l2/2; -h/2];
offset4 = [-l1/2; +l2/2; -h/2];
base1_0 = base_0 + R_SB*offset1;
base2_0 = base_0 + R_SB*offset2;
base3_0 = base_0 + R_SB*offset3;
base4_0 = base_0 + R_SB*offset4;

% First leg
theta_leg1_0 = [pi/6; pi/6; pi/6];
p_leg1_0 = forward_kin(theta_leg1_0);
T_SB1 = T_ypr(ypr_0(1), ypr_0(2), ypr_0(3), base1_0(1), base1_0(2), base1_0(3));
ee_leg1_0 = (T_SB1 * [p_leg1_0; 1]);
ee_leg1_0 = ee_leg1_0(1:3);

% Second leg
theta_leg2_0 = [pi/6; pi/6; pi/6];
p_leg2_0 = forward_kin(theta_leg2_0);
T_SB2 = T_ypr(ypr_0(1), ypr_0(2), ypr_0(3), base2_0(1), base2_0(2), base2_0(3));
ee_leg2_0 = double(T_SB2 * [p_leg2_0; 1]);
ee_leg2_0 = ee_leg2_0(1:3);

% Third leg
theta_leg3_0 = [pi/6; pi/6; pi/6];
p_leg3_0 = forward_kin(theta_leg3_0);
T_SB3 = T_ypr(ypr_0(1), ypr_0(2), ypr_0(3), base3_0(1), base3_0(2), base3_0(3));
ee_leg3_0 = double(T_SB3 * [p_leg3_0; 1]);
ee_leg3_0 = ee_leg3_0(1:3);

% Fourth leg
theta_leg4_0 = [pi/6; pi/6; pi/6];
p_leg4_0 = forward_kin(theta_leg4_0);
T_SB4 = T_ypr(ypr_0(1), ypr_0(2), ypr_0(3), base4_0(1), base4_0(2), base4_0(3));
ee_leg4_0 = double(T_SB4 * [p_leg4_0; 1]);
ee_leg4_0 = ee_leg4_0(1:3);

% Robot
q_0 = [base_0; ypr_0; theta_leg1_0; ee_leg1_0; theta_leg2_0; ee_leg2_0;...
       theta_leg3_0; ee_leg3_0; theta_leg4_0; ee_leg4_0;];

ni_0 = zeros(18,1);

%% Gait parameters: COG Constant Velocity
init_DST_const_vel;
ratio = 2;             % ratio shift_duration/swing_duration
H = 0.05;              % height of a swing
L = 0.3;               % length of a cycle of four swings
vel = 0.2;             % CoM velocity 
steps = 4+2*ratio;     % # of swings and shifts
T = L/vel;             % time to do a cycle of four swings
ts = T/steps;          % time to do a single swing 
l = L*(steps-1)/steps; % length of a single swing
ts2 = 3*ts;
% Points for rotation
points_cw = [-0.1,  0.00;
	     	+0.1, -0.05;
	      	0.0,     0];
points_ccw = [+0.1,  0.00;
	      	-0.1, -0.05;
	       	0.0,     0];
	
%% Guidance
init_guidance;

%% Controller 
K_p = 200;
K_d = 10; 
K_i = 1; 

%% Ground
K_el = 10000;
K_v = 10;
mu = 0.8;
