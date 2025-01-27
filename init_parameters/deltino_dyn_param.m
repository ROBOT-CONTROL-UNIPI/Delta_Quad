%% Robot dynamics parameters
m_t = 0.2;  % mass of the thigh 200 g
m_j = 0.1;  % mass of the knee joint 100 g
m_s = 0.15; % mass of the shin 150 g
m_f = 0.3;  % mass of the foot 300 g
m_p = 0.1;  % mass of the payload 0.1 g
m_b = 3;    % mass of the base 3.0 kg
m_c = m_f + m_p + 1/3*m_s;    % total mass of the coupler
m_k = m_j + 2/3*m_s;          % total mass of the knee
m_l = m_t + m_k;              % total mass of the upper leg
m_leg = 3*m_t + 3*m_k + m_c;  % total mass of one delta leg
m_tot = m_leg*4 + m_b;        % total mass of robot