function Tij = T_ij3D(x, y, n)
    r_vec = y-x;
    r = norm(r_vec);
    r_comma_1 = r_vec(1) / r;
    r_comma_2 = r_vec(2) / r;
    r_comma_3 = r_vec(3) / r;
    r_comma_k_n_k = r_comma_1*n(1) + r_comma_2*n(2) + r_comma_3*n(3);
    cst = -3/(4*pi * r^2) * r_comma_k_n_k;
    xyyx = r_comma_1 * r_comma_2;
    xzzx = r_comma_1 * r_comma_3;
    zyyz = r_comma_2 * r_comma_3;
    Tij = cst*[r_comma_1^2, xyyx, xzzx; xyyx, r_comma_2^2, zyyz; xzzx, zyyz, r_comma_3^2];
end
