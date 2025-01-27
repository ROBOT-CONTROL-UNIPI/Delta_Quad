%% Robot geometric parameters
l_A = 0.2;      % thigh length 
l_B = 0.4;      % shin length 
l1 = 0.7;       % length of central body
l2 = 0.7;       % width of central body 
h = 0.2;        % height of central body
R_A = 0.15;     % base plate radius 
R_B = 0.05;     % travelling plate radius 
R = R_A - R_B;

% Orientation of single frames {B_i} w.r.t body frame {B}
alpha1 = 0;
alpha2 = 2/3*pi;
alpha3 = 4/3*pi;

alpha11 = +pi/4;
alpha12 = 2/3*pi+pi/4;
alpha13 = 4/3*pi+pi/4;

alpha21 = -pi/4;
alpha22 = 2/3*pi-pi/4;
alpha23 = 4/3*pi-pi/4; 

alpha31 = -pi/12;
alpha32 = 2/3*pi-pi/12;
alpha33 = 4/3*pi-pi/12; 

alpha41 = +pi/12;
alpha42 = 2/3*pi+pi/12;
alpha43 = 4/3*pi+pi/12;