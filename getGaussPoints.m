function [gp_x, gw_x, gp_y, gw_y] = getGaussPoints(ngp_x, ngp_y)
    % For 3D surface integral we need 2 gauss pairs(eta, epsilon) and 2 integrals.(Heinz Antes a short course on boundary element methods book)
    [gw_x, gp_x] = gauss(ngp_x);
    [gw_y, gp_y] = gauss(ngp_y);
end
