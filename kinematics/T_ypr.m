function T = T_ypr(phi, theta, psi, x, y, z)

T = [cos(phi)*cos(theta), cos(phi)*sin(psi)*sin(theta) - cos(psi)*sin(phi), sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta), x;
     cos(theta)*sin(phi), cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta), cos(psi)*sin(phi)*sin(theta) - cos(phi)*sin(psi), y;
             -sin(theta),                              cos(theta)*sin(psi),                              cos(psi)*cos(theta), z;
                       0,                                                0,                                                0, 1];
end