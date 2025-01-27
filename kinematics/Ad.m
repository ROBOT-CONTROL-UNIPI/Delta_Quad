function Ad_g = Ad(g)

    R = g(1:3, 1:3);
    d = g(1:3, 4);

    Ad_g = [       R, hat(d)*R;
            zeros(3),        R];

end