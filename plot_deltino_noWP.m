%% Init
clc; close gcf; 

scale = 1; % scale factor for plotting rate
dt = 1/60; % sampling time for plotting [s]

deltino_kin_param

%% Obtaining Simulink data
model = 'deltino_feet_ref_stateflow';

% fprintf('Loading Simulink model...\n')
% load_system(model);
% fprintf('Running Simulink simulation ...\n')
% out = sim(model);

% Reorganizing Simulink timeseries
fprintf('Obtaining data from Simulink ...\n')
dim = length(out.q.Time);
tout = out.tout;
base = reshape(out.q.Data(1:3, :, :), 3, dim);
ypr = reshape(out.q.Data(4:6, :, :), 3, dim);
theta1 = reshape(out.q.Data(7:9, :, :), 3, dim);
ee1 = reshape(out.q.Data(10:12, :, :), 3, dim);
theta2 = reshape(out.q.Data(13:15, :, :), 3, dim);
ee2 = reshape(out.q.Data(16:18, :, :), 3, dim);
theta3 = reshape(out.q.Data(19:21, :, :), 3, dim);
ee3 = reshape(out.q.Data(22:24, :, :), 3, dim);
theta4 = reshape(out.q.Data(25:27, :, :), 3, dim);
ee4 = reshape(out.q.Data(28:30, :, :), 3, dim);
fprintf('... data from Simulink obtained\n')

% Post processing data
fprintf('Post processing of data ...\n')
time = 0 : dt : tout(end);
steps2 = length(time);
base_pp = interp1(tout, base.', time).';
ypr_pp = interp1(tout, ypr.', time).';
theta1_pp = interp1(tout, theta1.', time).';
ee1_pp = interp1(tout, ee1.', time).';
theta2_pp = interp1(tout, theta2.', time).';
ee2_pp = interp1(tout, ee2.', time).';
theta3_pp = interp1(tout, theta3.', time).';
ee3_pp = interp1(tout, ee3.', time).';
theta4_pp = interp1(tout, theta4.', time).';
ee4_pp = interp1(tout, ee4.', time).';
fprintf('... post processing terminated\n')

%% Computing base, hips, knees and ankles coordinates of first delta in frame {S}
% Hip joints coordinates in frame {S}
base1 = zeros(3, steps2);
Pa11 = zeros(4,steps2);
Pa12 = zeros(4,steps2);
Pa13 = zeros(4,steps2);
for i=1:steps2
    R_SB = Rypr(ypr_pp(1,i), ypr_pp(2,i), ypr_pp(3,i));
    base1(:,i) = base_pp(:,i) + R_SB * [l1/2; l2/2;-h/2];  
    T_SB1 = T_ypr(ypr_pp(1,i), ypr_pp(2,i), ypr_pp(3,i), base1(1,i), base1(2,i), base1(3,i));
    Pa11(:,i) = T_SB1 * R_z(alpha1) * [R_A, 0, 0, 1].';
    Pa12(:,i) = T_SB1 * R_z(alpha2) * [R_A, 0, 0, 1].';
    Pa13(:,i) = T_SB1 * R_z(alpha3) * [R_A, 0, 0, 1].';
end 
% Knees coordinates in frame {S}
Pc11 = zeros(4,steps2);
Pc12 = zeros(4,steps2);
Pc13 = zeros(4,steps2);
for i = 1:steps2
    T_SB1 = T_ypr(ypr_pp(1,i), ypr_pp(2,i), ypr_pp(3,i), base1(1,i), base1(2,i), base1(3,i));
    Pc11(:,i) = T_SB1 * R_z(alpha1) * [R_A+l_A*cos(theta1_pp(1,i)), 0, -l_A*sin(theta1_pp(1,i)), 1].';
    Pc12(:,i) = T_SB1 * R_z(alpha2) * [R_A+l_A*cos(theta1_pp(2,i)), 0, -l_A*sin(theta1_pp(2,i)), 1].';
    Pc13(:,i) = T_SB1 * R_z(alpha3) * [R_A+l_A*cos(theta1_pp(3,i)), 0, -l_A*sin(theta1_pp(3,i)), 1].';
end
% Foot coordinates in frame {S}
P1 = zeros(4,steps2);
for i=1:steps2
    P1(:,i) = [ee1_pp(:,i); 1];
end
% Ankles coordinates in frame {S}
Pb11 = zeros(4,steps2);
Pb12 = zeros(4,steps2);
Pb13 = zeros(4,steps2);
for i = 1:steps2
    R_SB = T_ypr(ypr_pp(1,i), ypr_pp(2,i), ypr_pp(3,i), 0, 0, 0);
    Pb11(:,i) = P1(:,i) + R_SB * [R_B*cos(alpha1), R_B*sin(alpha1), 0, 1].' - [0;0;0;1];
    Pb12(:,i) = P1(:,i) + R_SB * [R_B*cos(alpha2), R_B*sin(alpha2), 0, 1].' - [0;0;0;1];
    Pb13(:,i) = P1(:,i) + R_SB * [R_B*cos(alpha3), R_B*sin(alpha3), 0, 1].' - [0;0;0;1];
end

fprintf('Computing of first leg coordinates terminated\n')

%% Computing base, hips, knees and ankles coordinates of second delta in frame {S}
% Hip joints coordinates in frame {S}
base2 = zeros(3, steps2);
Pa21 = zeros(4,steps2);
Pa22 = zeros(4,steps2);
Pa23 = zeros(4,steps2);
for i=1:steps2
    R_SB = Rypr(ypr_pp(1,i), ypr_pp(2,i), ypr_pp(3,i));
    base2(:,i) = base_pp(:,i) + R_SB * [l1/2; -l2/2;-h/2];  
    T_SB2 = T_ypr(ypr_pp(1,i), ypr_pp(2,i), ypr_pp(3,i), base2(1,i), base2(2,i), base2(3,i));
    Pa21(:,i) = T_SB2 * R_z(alpha1) * [R_A, 0, 0, 1].';
    Pa22(:,i) = T_SB2 * R_z(alpha2) * [R_A, 0, 0, 1].';
    Pa23(:,i) = T_SB2 * R_z(alpha3) * [R_A, 0, 0, 1].';
end 
% Knees coordinates in frame {S}
Pc21 = zeros(4,steps2);
Pc22 = zeros(4,steps2);
Pc23 = zeros(4,steps2);
for i = 1:steps2
    T_SB2 = T_ypr(ypr_pp(1,i), ypr_pp(2,i), ypr_pp(3,i), base2(1,i), base2(2,i), base2(3,i));
    Pc21(:,i) = T_SB2 * R_z(alpha1) * [R_A+l_A*cos(theta2_pp(1,i)), 0, -l_A*sin(theta2_pp(1,i)), 1].';
    Pc22(:,i) = T_SB2 * R_z(alpha2) * [R_A+l_A*cos(theta2_pp(2,i)), 0, -l_A*sin(theta2_pp(2,i)), 1].';
    Pc23(:,i) = T_SB2 * R_z(alpha3) * [R_A+l_A*cos(theta2_pp(3,i)), 0, -l_A*sin(theta2_pp(3,i)), 1].';
end
% Foot coordinates in frame {S}
P2 = zeros(4,steps2);
for i=1:steps2
    P2(:,i) = [ee2_pp(:,i); 1];
end
% Ankles coordinates in frame {S}
Pb21 = zeros(4,steps2);
Pb22 = zeros(4,steps2);
Pb23 = zeros(4,steps2);
for i = 1:steps2
    R_SB = T_ypr(ypr_pp(1,i), ypr_pp(2,i), ypr_pp(3,i), 0, 0, 0);
    Pb21(:,i) = P2(:,i) + R_SB * [R_B*cos(alpha1), R_B*sin(alpha1), 0, 1].' - [0;0;0;1];
    Pb22(:,i) = P2(:,i) + R_SB * [R_B*cos(alpha2), R_B*sin(alpha2), 0, 1].' - [0;0;0;1];
    Pb23(:,i) = P2(:,i) + R_SB * [R_B*cos(alpha3), R_B*sin(alpha3), 0, 1].' - [0;0;0;1];
end

fprintf('Computing of second leg coordinates terminated\n')

%% Computing base, hips, knees and ankles coordinates of third delta in frame {S}
% Hip joints coordinates in frame {S}
base3 = zeros(3, steps2);
Pa31 = zeros(4,steps2);
Pa32 = zeros(4,steps2);
Pa33 = zeros(4,steps2);
for i=1:steps2
    R_SB = Rypr(ypr_pp(1,i), ypr_pp(2,i), ypr_pp(3,i));
    base3(:,i) = base_pp(:,i) + R_SB * [-l1/2; -l2/2;-h/2];  
    T_SB3 = T_ypr(ypr_pp(1,i), ypr_pp(2,i), ypr_pp(3,i), base3(1,i), base3(2,i), base3(3,i));
    Pa31(:,i) = T_SB3 * R_z(alpha1) * [R_A, 0, 0, 1].';
    Pa32(:,i) = T_SB3 * R_z(alpha2) * [R_A, 0, 0, 1].';
    Pa33(:,i) = T_SB3 * R_z(alpha3) * [R_A, 0, 0, 1].';
end 
% Knees coordinates in frame {S}
Pc31 = zeros(4,steps2);
Pc32 = zeros(4,steps2);
Pc33 = zeros(4,steps2);
for i = 1:steps2
    T_SB3 = T_ypr(ypr_pp(1,i), ypr_pp(2,i), ypr_pp(3,i), base3(1,i), base3(2,i), base3(3,i));
    Pc31(:,i) = T_SB3 * R_z(alpha1) * [R_A+l_A*cos(theta3_pp(1,i)), 0, -l_A*sin(theta3_pp(1,i)), 1].';
    Pc32(:,i) = T_SB3 * R_z(alpha2) * [R_A+l_A*cos(theta3_pp(2,i)), 0, -l_A*sin(theta3_pp(2,i)), 1].';
    Pc33(:,i) = T_SB3 * R_z(alpha3) * [R_A+l_A*cos(theta3_pp(3,i)), 0, -l_A*sin(theta3_pp(3,i)), 1].';
end
% Foot coordinates in frame {S}
P3 = zeros(4,steps2);
for i=1:steps2
    P3(:,i) = [ee3_pp(:,i); 1];
end
% Ankles coordinates in frame {S}
Pb31 = zeros(4,steps2);
Pb32 = zeros(4,steps2);
Pb33 = zeros(4,steps2);
for i = 1:steps2
    R_SB = T_ypr(ypr_pp(1,i), ypr_pp(2,i), ypr_pp(3,i), 0, 0, 0);
    Pb31(:,i) = P3(:,i) + R_SB * [R_B*cos(alpha1), R_B*sin(alpha1), 0, 1].' - [0;0;0;1];
    Pb32(:,i) = P3(:,i) + R_SB * [R_B*cos(alpha2), R_B*sin(alpha2), 0, 1].' - [0;0;0;1];
    Pb33(:,i) = P3(:,i) + R_SB * [R_B*cos(alpha3), R_B*sin(alpha3), 0, 1].' - [0;0;0;1];
end

fprintf('Computing of third leg coordinates terminated\n')

%% Computing base, hips, knees and ankles coordinates of forth delta in frame {S}
% Hip joints coordinates in frame {S}
base4 = zeros(3, steps2);
Pa41 = zeros(4,steps2);
Pa42 = zeros(4,steps2);
Pa43 = zeros(4,steps2);
for i=1:steps2
    R_SB = Rypr(ypr_pp(1,i), ypr_pp(2,i), ypr_pp(3,i));
    base4(:,i) = base_pp(:,i) + R_SB * [-l1/2; l2/2; -h/2];  
    T_SB4 = T_ypr(ypr_pp(1,i), ypr_pp(2,i), ypr_pp(3,i), base4(1,i), base4(2,i), base4(3,i));
    Pa41(:,i) = T_SB4 * R_z(alpha1) * [R_A, 0, 0, 1].';
    Pa42(:,i) = T_SB4 * R_z(alpha2) * [R_A, 0, 0, 1].';
    Pa43(:,i) = T_SB4 * R_z(alpha3) * [R_A, 0, 0, 1].';
end 
% Knees coordinates in frame {S}
Pc41 = zeros(4,steps2);
Pc42 = zeros(4,steps2);
Pc43 = zeros(4,steps2);
for i = 1:steps2
    T_SB4 = T_ypr(ypr_pp(1,i), ypr_pp(2,i), ypr_pp(3,i), base4(1,i), base4(2,i), base4(3,i));
    Pc41(:,i) = T_SB4 * R_z(alpha1) * [R_A+l_A*cos(theta4_pp(1,i)), 0, -l_A*sin(theta4_pp(1,i)), 1].';
    Pc42(:,i) = T_SB4 * R_z(alpha2) * [R_A+l_A*cos(theta4_pp(2,i)), 0, -l_A*sin(theta4_pp(2,i)), 1].';
    Pc43(:,i) = T_SB4 * R_z(alpha3) * [R_A+l_A*cos(theta4_pp(3,i)), 0, -l_A*sin(theta4_pp(3,i)), 1].';
end
% Foot coordinates in frame {S}
P4 = zeros(4,steps2);
for i=1:steps2
    P4(:,i) = [ee4_pp(:,i); 1];
end
% Ankles coordinates in frame {S}
Pb41 = zeros(4,steps2);
Pb42 = zeros(4,steps2);
Pb43 = zeros(4,steps2);
for i = 1:steps2
    R_SB = T_ypr(ypr_pp(1,i), ypr_pp(2,i), ypr_pp(3,i), 0, 0, 0);
    Pb41(:,i) = P4(:,i) + R_SB * [R_B*cos(alpha1), R_B*sin(alpha1), 0, 1].' - [0;0;0;1];
    Pb42(:,i) = P4(:,i) + R_SB * [R_B*cos(alpha2), R_B*sin(alpha2), 0, 1].' - [0;0;0;1];
    Pb43(:,i) = P4(:,i) + R_SB * [R_B*cos(alpha3), R_B*sin(alpha3), 0, 1].' - [0;0;0;1];
end

fprintf('Computing of forth leg coordinates terminated\n')

%% Initial plot
% Plot spatial frame
plot3([0, 0.2], [0, 0], [0, 0], 'Color', 'k', 'Linewidth', 1);
hold on
text(0.2, 0, 0.05, "X");
plot3([0, 0], [0, 0.2], [0, 0], 'Color', 'k', 'Linewidth', 1);
text(0, 0.2, 0.05, "Y");
plot3([0, 0], [0, 0], [0, 0.2], 'Color', 'k', 'Linewidth', 1);
text(0, 0, 0.25, "Z");

% Other graphical features
grid on
axis equal
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
xlim([-1.5, 1.7]);
ylim([-1.5, 1.5]); 
zlim([-0.05, 1]);
view([-90,0])

% % Floor plot
% [x, y] = meshgrid(-10:0.1:10); % Generate x and y data
% z = zeros(size(x, 1)); % Generate z data
% surf(x, y, z, 'EdgeColor', 'none', 'FaceAlpha', '0.0')

F = struct('cdata', cell(1,steps2), 'colormap', cell(1,steps2));
duration = zeros(steps2,1);
fig = gcf;
fig.WindowState = 'maximized';

%% Cyclic plot of robot
for i = 1:steps2
    tic
   
    % Plot support triangle
    idx = 0;
    support_x = zeros(5,1);
    support_y = zeros(5,1);
    if P1(3,i)<=0.005
        idx = idx +1;
        support_x(idx) = P1(1,i);
        support_y(idx) = P1(2,i);
    end
    if P2(3,i)<=0.005
        idx = idx +1;
        support_x(idx) = P2(1,i);
        support_y(idx) = P2(2,i);
    end
    if P3(3,i)<=0.005
        idx = idx +1;
        support_x(idx) = P3(1,i);
        support_y(idx) = P3(2,i);
    end
    if P4(3,i)<=0.005
        idx = idx +1;
        support_x(idx) = P4(1,i);
        support_y(idx) = P4(2,i);
	end
    idx = idx +1;
    support_x(idx) = support_x(1);
    support_y(idx) = support_y(1);
    support_x = support_x(1:idx);
    support_y = support_y(1:idx);
    res0_2 = fill3(support_x, support_y, zeros(1,idx), 'y');
    res0_2.FaceAlpha = 0.5; 
    % Plot CoM projection
    CoM_1 = (m_t/2*Pa11(1:3,i) + m_t/2*Pa12(1:3,i) + m_t/2*Pa13(1:3,i) + ... % hips
            (m_t/2+m_k)*Pc11(1:3,i) + (m_t/2+m_k)*Pc12(1:3,i) + (m_t/2+m_k)*Pc13(1:3,i) + ... % knees
            m_c*P1(1:3,i))/m_leg;% coupler
    CoM_2 = (m_t/2*Pa21(1:3,i) + m_t/2*Pa22(1:3,i) + m_t/2*Pa23(1:3,i) + ... % hips
            (m_t/2+m_k)*Pc21(1:3,i) + (m_t/2+m_k)*Pc22(1:3,i) + (m_t/2+m_k)*Pc23(1:3,i) + ... % knees
            m_c*P2(1:3,i))/m_leg;% coupler
    CoM_3 = (m_t/2*Pa31(1:3,i) + m_t/2*Pa32(1:3,i) + m_t/2*Pa33(1:3,i) + ... % hips
            (m_t/2+m_k)*Pc31(1:3,i) + (m_t/2+m_k)*Pc32(1:3,i) + (m_t/2+m_k)*Pc33(1:3,i) + ... % knees
            m_c*P3(1:3,i))/m_leg;% coupler
    CoM_4 = (m_t/2*Pa41(1:3,i) + m_t/2*Pa42(1:3,i) + m_t/2*Pa43(1:3,i) + ... % hips
            (m_t/2+m_k)*Pc41(1:3,i) + (m_t/2+m_k)*Pc42(1:3,i) + (m_t/2+m_k)*Pc43(1:3,i) + ... % knees
            m_c*P4(1:3,i))/m_leg;% coupler
    CoM = (CoM_1*m_leg + CoM_2*m_leg + CoM_3*m_leg + CoM_4*m_leg + base_pp(1:3,i)*m_b)/m_tot;
    res0_3 = plot3(base_pp(1,i), base_pp(2,i), 0, '*', 'Color', 'red');
    res0_4 = plot3(CoM(1), CoM(2), 0, '*', 'Color', 'green');

    % Plot base of first delta
    res1_1 = fill3([Pa11(1,i); Pa12(1,i); Pa13(1,i); Pa11(1,i)], ...
                   [Pa11(2,i); Pa12(2,i); Pa13(2,i); Pa11(2,i)], ...
                   [Pa11(3,i); Pa12(3,i); Pa13(3,i); Pa11(3,i)], ...
                   [0.5 0.5 0.5]);
	res1_1.LineWidth = 1;
	res1_1.FaceAlpha = 0.95;
    % Plot hip joints of first delta
    res1_2 = plot3(Pa11(1,i),Pa11(2,i),Pa11(3,i), 'o', 'Color', 'k', 'LineWidth', 1);
    res1_3 = plot3(Pa12(1,i),Pa12(2,i),Pa12(3,i), 'o', 'Color', 'k', 'LineWidth', 1);
    res1_4 = plot3(Pa13(1,i),Pa13(2,i),Pa13(3,i), 'o', 'Color', 'k', 'LineWidth', 1);
    % Plot coupler of first delta
    res1_5 = plot3([Pb11(1,i), Pb12(1,i), Pb13(1,i), Pb11(1,i)], ...
                   [Pb11(2,i), Pb12(2,i), Pb13(2,i), Pb11(2,i)], ...
                   [Pb11(3,i), Pb12(3,i), Pb13(3,i), Pb11(3,i)], ...
                   'Color', 'k');
    % Plot legs of first delta
    res1_6 = plot3([P1(1,i), Pb11(1,i), Pc11(1,i), Pa11(1,i)], ...
                   [P1(2,i), Pb11(2,i), Pc11(2,i), Pa11(2,i)], ...
                   [P1(3,i), Pb11(3,i), Pc11(3,i), Pa11(3,i)], ...
                   'Color', 'b', 'LineWidth', 1);
    res1_7 = plot3([P1(1,i), Pb12(1,i), Pc12(1,i), Pa12(1,i)], ...
                   [P1(2,i), Pb12(2,i), Pc12(2,i), Pa12(2,i)], ...
                   [P1(3,i), Pb12(3,i), Pc12(3,i), Pa12(3,i)], ...
                   'Color', 'b', 'LineWidth', 1);
    res1_8 = plot3([P1(1,i), Pb13(1,i), Pc13(1,i), Pa13(1,i)], ...
                   [P1(2,i), Pb13(2,i), Pc13(2,i), Pa13(2,i)], ...
                   [P1(3,i), Pb13(3,i), Pc13(3,i), Pa13(3,i)], ...
                   'Color', 'b', 'LineWidth', 1);
    % Plot knee joints of first delta
    res1_9  = plot3(Pc11(1,i),Pc11(2,i),Pc11(3,i), 'o', 'Color', 'k', 'LineWidth', 1);
    res1_10 = plot3(Pc12(1,i),Pc12(2,i),Pc12(3,i), 'o', 'Color', 'k', 'LineWidth', 1);
    res1_11 = plot3(Pc13(1,i),Pc13(2,i),Pc13(3,i), 'o', 'Color', 'k', 'LineWidth', 1);
    % Plot ankle joints of first delta
    res1_12 = plot3(Pb11(1,i),Pb11(2,i),Pb11(3,i), 'o', 'Color', 'k', 'LineWidth', 1);
    res1_13 = plot3(Pb12(1,i),Pb12(2,i),Pb12(3,i), 'o', 'Color', 'k', 'LineWidth', 1);
    res1_14 = plot3(Pb13(1,i),Pb13(2,i),Pb13(3,i), 'o', 'Color', 'k', 'LineWidth', 1);

    % Plot base of second delta
	res2_1 = fill3([Pa21(1,i); Pa22(1,i); Pa23(1,i); Pa21(1,i)], ...
                   [Pa21(2,i); Pa22(2,i); Pa23(2,i); Pa21(2,i)], ...
                   [Pa21(3,i); Pa22(3,i); Pa23(3,i); Pa21(3,i)], ...
                   [0.5 0.5 0.5]);
	res2_1.LineWidth = 1;
	res2_1.FaceAlpha = 0.95;
    % Plot hip joints of second delta
    res2_2 = plot3(Pa21(1,i),Pa21(2,i),Pa21(3,i), 'o', 'Color', 'k', 'LineWidth', 1);
    res2_3 = plot3(Pa22(1,i),Pa22(2,i),Pa22(3,i), 'o', 'Color', 'k', 'LineWidth', 1);
    res2_4 = plot3(Pa23(1,i),Pa23(2,i),Pa23(3,i), 'o', 'Color', 'k', 'LineWidth', 1);
    % Plot coupler of second delta
    res2_5 = plot3([Pb21(1,i), Pb22(1,i), Pb23(1,i), Pb21(1,i)], ...
                   [Pb21(2,i), Pb22(2,i), Pb23(2,i), Pb21(2,i)], ...
                   [Pb21(3,i), Pb22(3,i), Pb23(3,i), Pb21(3,i)], ...
                   'Color', 'k', 'LineWidth', 1);
    % Plot legs of second delta
    res2_6 = plot3([P2(1,i), Pb21(1,i), Pc21(1,i), Pa21(1,i)], ...
                   [P2(2,i), Pb21(2,i), Pc21(2,i), Pa21(2,i)], ...
                   [P2(3,i), Pb21(3,i), Pc21(3,i), Pa21(3,i)], ...
                   'Color', 'g', 'LineWidth', 1);
    res2_7 = plot3([P2(1,i), Pb22(1,i), Pc22(1,i), Pa22(1,i)], ...
                   [P2(2,i), Pb22(2,i), Pc22(2,i), Pa22(2,i)], ...
                   [P2(3,i), Pb22(3,i), Pc22(3,i), Pa22(3,i)], ...
                   'Color', 'g', 'LineWidth', 1);
    res2_8 = plot3([P2(1,i), Pb23(1,i), Pc23(1,i), Pa23(1,i)], ...
                   [P2(2,i), Pb23(2,i), Pc23(2,i), Pa23(2,i)], ...
                   [P2(3,i), Pb23(3,i), Pc23(3,i), Pa23(3,i)], ...
                   'Color', 'g', 'LineWidth', 1);
    % Plot knee joints of second delta
    res2_9  = plot3(Pc21(1,i),Pc21(2,i),Pc21(3,i), 'o', 'Color', 'k', 'LineWidth', 1);
    res2_10 = plot3(Pc22(1,i),Pc22(2,i),Pc22(3,i), 'o', 'Color', 'k', 'LineWidth', 1);
    res2_11 = plot3(Pc23(1,i),Pc23(2,i),Pc23(3,i), 'o', 'Color', 'k', 'LineWidth', 1);
    % Plot ankle joints of second delta
    res2_12 = plot3(Pb21(1,i),Pb21(2,i),Pb21(3,i), 'o', 'Color', 'k', 'LineWidth', 1);
    res2_13 = plot3(Pb22(1,i),Pb22(2,i),Pb22(3,i), 'o', 'Color', 'k', 'LineWidth', 1);
    res2_14 = plot3(Pb23(1,i),Pb23(2,i),Pb23(3,i), 'o', 'Color', 'k', 'LineWidth', 1); 

    % Plot base of third delta
    res3_1 = fill3([Pa31(1,i); Pa32(1,i); Pa33(1,i); Pa31(1,i)], ...
                   [Pa31(2,i); Pa32(2,i); Pa33(2,i); Pa31(2,i)], ...
                   [Pa31(3,i); Pa32(3,i); Pa33(3,i); Pa31(3,i)], ...
                   [0.5 0.5 0.5]);
	res3_1.LineWidth = 1;
	res3_1.FaceAlpha = 0.95;
    % Plot hip joints of third delta
    res3_2 = plot3(Pa31(1,i),Pa31(2,i),Pa31(3,i), 'o', 'Color', 'k', 'LineWidth', 1);
    res3_3 = plot3(Pa32(1,i),Pa32(2,i),Pa32(3,i), 'o', 'Color', 'k', 'LineWidth', 1);
    res3_4 = plot3(Pa33(1,i),Pa33(2,i),Pa33(3,i), 'o', 'Color', 'r', 'LineWidth', 1);
    % Plot coupler of third delta
    res3_5 = plot3([Pb31(1,i), Pb32(1,i), Pb33(1,i), Pb31(1,i)], ...
                   [Pb31(2,i), Pb32(2,i), Pb33(2,i), Pb31(2,i)], ...
                   [Pb31(3,i), Pb32(3,i), Pb33(3,i), Pb31(3,i)], ...
                   'Color', 'k');
    % Plot legs of third delta
    res3_6 = plot3([P3(1,i), Pb31(1,i), Pc31(1,i), Pa31(1,i)], ...
                   [P3(2,i), Pb31(2,i), Pc31(2,i), Pa31(2,i)], ...
                   [P3(3,i), Pb31(3,i), Pc31(3,i), Pa31(3,i)], ...
                   'Color', 'r', 'LineWidth', 1);
    res3_7 = plot3([P3(1,i), Pb32(1,i), Pc32(1,i), Pa32(1,i)], ...
                   [P3(2,i), Pb32(2,i), Pc32(2,i), Pa32(2,i)], ...
                   [P3(3,i), Pb32(3,i), Pc32(3,i), Pa32(3,i)], ...
                   'Color', 'r', 'LineWidth', 1);
    res3_8 = plot3([P3(1,i), Pb33(1,i), Pc33(1,i), Pa33(1,i)], ...
                   [P3(2,i), Pb33(2,i), Pc33(2,i), Pa33(2,i)], ...
                   [P3(3,i), Pb33(3,i), Pc33(3,i), Pa33(3,i)], ...
                   'Color', 'r', 'LineWidth', 1);
    % Plot knee joints of third delta
    res3_9  = plot3(Pc31(1,i),Pc31(2,i),Pc31(3,i), 'o', 'Color', 'k', 'LineWidth', 1);
    res3_10 = plot3(Pc32(1,i),Pc32(2,i),Pc32(3,i), 'o', 'Color', 'k', 'LineWidth', 1);
    res3_11 = plot3(Pc33(1,i),Pc33(2,i),Pc33(3,i), 'o', 'Color', 'k', 'LineWidth', 1);
    % Plot ankle joints of third delta
    res3_12 = plot3(Pb31(1,i),Pb31(2,i),Pb31(3,i), 'o', 'Color', 'k', 'LineWidth', 1);
    res3_13 = plot3(Pb32(1,i),Pb32(2,i),Pb32(3,i), 'o', 'Color', 'k', 'LineWidth', 1);
    res3_14 = plot3(Pb33(1,i),Pb33(2,i),Pb33(3,i), 'o', 'Color', 'k', 'LineWidth', 1);
 
    % Plot base of forth delta
	res4_1 = fill3([Pa41(1,i); Pa42(1,i); Pa43(1,i); Pa41(1,i)], ...
                   [Pa41(2,i); Pa42(2,i); Pa43(2,i); Pa41(2,i)], ...
                   [Pa41(3,i); Pa42(3,i); Pa43(3,i); Pa41(3,i)], ...
                   [0.5 0.5 0.5]);
	res4_1.LineWidth = 1;
	res4_1.FaceAlpha = 0.95;
    % Plot hip joints of forth delta
    res4_2 = plot3(Pa41(1,i),Pa41(2,i),Pa41(3,i), 'o', 'Color', 'k', 'LineWidth', 1);
    res4_3 = plot3(Pa42(1,i),Pa42(2,i),Pa42(3,i), 'o', 'Color', 'k', 'LineWidth', 1);
    res4_4 = plot3(Pa43(1,i),Pa43(2,i),Pa43(3,i), 'o', 'Color', 'k', 'LineWidth', 1);
    % Plot coupler of forth delta
    res4_5 = plot3([Pb41(1,i), Pb42(1,i), Pb43(1,i), Pb41(1,i)], ...
                   [Pb41(2,i), Pb42(2,i), Pb43(2,i), Pb41(2,i)], ...
                   [Pb41(3,i), Pb42(3,i), Pb43(3,i), Pb41(3,i)], ...
                   'Color', 'k', 'LineWidth', 1);
    % Plot legs of forth delta
    res4_6 = plot3([P4(1,i), Pb41(1,i), Pc41(1,i), Pa41(1,i)], ...
                   [P4(2,i), Pb41(2,i), Pc41(2,i), Pa41(2,i)], ...
                   [P4(3,i), Pb41(3,i), Pc41(3,i), Pa41(3,i)], ...
                   'Color', 'cyan', 'LineWidth', 1);
    res4_7 = plot3([P4(1,i), Pb42(1,i), Pc42(1,i), Pa42(1,i)], ...
                   [P4(2,i), Pb42(2,i), Pc42(2,i), Pa42(2,i)], ...
                   [P4(3,i), Pb42(3,i), Pc42(3,i), Pa42(3,i)], ...
                   'Color', 'cyan', 'LineWidth', 1);
    res4_8 = plot3([P4(1,i), Pb43(1,i), Pc43(1,i), Pa43(1,i)], ...
                   [P4(2,i), Pb43(2,i), Pc43(2,i), Pa43(2,i)], ...
                   [P4(3,i), Pb43(3,i), Pc43(3,i), Pa43(3,i)], ...
                   'Color', 'cyan', 'LineWidth', 1);
    % Plot knee joints of forth delta
    res4_9  = plot3(Pc41(1,i),Pc41(2,i),Pc41(3,i), 'o', 'Color', 'k', 'LineWidth', 1);
    res4_10 = plot3(Pc42(1,i),Pc42(2,i),Pc42(3,i), 'o', 'Color', 'k', 'LineWidth', 1);
    res4_11 = plot3(Pc43(1,i),Pc43(2,i),Pc43(3,i), 'o', 'Color', 'k', 'LineWidth', 1);
    % Plot ankle joints of forth delta
    res4_12 = plot3(Pb41(1,i),Pb41(2,i),Pb41(3,i), 'o', 'Color', 'k', 'LineWidth', 1);
    res4_13 = plot3(Pb42(1,i),Pb42(2,i),Pb42(3,i), 'o', 'Color', 'k', 'LineWidth', 1);
    res4_14 = plot3(Pb43(1,i),Pb43(2,i),Pb43(3,i), 'o', 'Color', 'k', 'LineWidth', 1);

	% Plot simulation clock
    res0 = text(0, -0.05, 0, ["Time:" num2str(time(i))]);
    % Plot trunk of robot
    res0_1 = fill3([base1(1,i);base2(1,i);base3(1,i);base4(1,i);base1(1,i)], ...
                   [base1(2,i);base2(2,i);base3(2,i);base4(2,i);base1(2,i)], ...
                   [base1(3,i);base2(3,i);base3(3,i);base4(3,i);base1(3,i)],[0.5 0.5 0.5]);
	res0_1.FaceAlpha = 0.95; 
	res0_1.LineWidth = 1;

    if i == steps2
        break;
    end

    res = [res0, res0_1, res0_2, res0_3, res0_4, ... 
           res1_1, res1_2, res1_3, res1_4, res1_5, res1_6, res1_7, res1_8, res1_9, res1_10, res1_11, res1_12, res1_13, res1_14, ...
           res2_1, res2_2, res2_3, res2_4, res2_5, res2_6, res2_7, res2_8, res2_9, res2_10, res2_11, res2_12, res2_13, res2_14, ...
           res3_1, res3_2, res3_3, res3_4, res3_5, res3_6, res3_7, res3_8, res3_9, res3_10, res3_11, res3_12, res3_13, res3_14, ...
           res4_1, res4_2, res4_3, res4_4, res4_5, res4_6, res4_7, res4_8, res4_9, res4_10, res4_11, res4_12, res4_13, res4_14];
    
    %     F(i) = getframe(fig, [0 0 fig.Position(3) fig.Position(4)]);
    F(i) = getframe(fig, [fig.Position(1) fig.Position(2) ...
                          fig.Position(3) fig.Position(4)]);
    duration(i) = toc;
    if duration < dt*scale
        pause(dt*scale-duration)
    end

    delete(res);
end
close gcf
fig2 = figure;
fig2.WindowState = 'maximized';
pause();
movie(fig2, F, 1, 1/dt, [0 0 0 0]);