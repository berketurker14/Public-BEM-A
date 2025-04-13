function Uij = U_ij3D(x,y,viscosity)
    r_vec = y-x;
    r = norm(r_vec);
    r_comma_1 = r_vec(1)/r;
    r_comma_2 = r_vec(2)/r;
    r_comma_3 = r_vec(3)/r;
    multiplier = 1/(8*pi*viscosity*r);
    xyyx = r_comma_1 * r_comma_2;
    xzzx = r_comma_1 * r_comma_3;
    zyyz = r_comma_2 * r_comma_3;

    Uij = multiplier*[1+r_comma_1^2, xyyx, xzzx; xyyx, 1+r_comma_2^2, zyyz; xzzx, zyyz, 1+r_comma_3^2];
end
