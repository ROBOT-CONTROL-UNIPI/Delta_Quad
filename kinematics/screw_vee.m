function v = screw_vee(M)

v = [vee(M(1:3, 1:3));
            M(1:3, 4)];
end