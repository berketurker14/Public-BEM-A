function [normv_all, x_gauss_all] = precomputeElementData(p, t, gp_x, gp_y, Lx, Ly, Lz) % Add detJ calculation here

totalElements = size(t,1);
ngp_x = length(gp_x);
ngp_y = length(gp_y);
ngp_total = ngp_x*ngp_y;

normv_all = zeros(totalElements, 3);
x_gauss_all = zeros(totalElements, ngp_total*3);

for elem = 1:totalElements
    elemNodes = t(elem,:);
    x1 = p(elemNodes(1),:);
    x2 = p(elemNodes(2),:);
    x3 = p(elemNodes(3),:);
    x4 = p(elemNodes(4),:);

    vector1 = x2 - x1;
    vector2 = x4 - x1;
    normal = cross(vector1, vector2);
    normv = normal / norm(normal);

    desiredDirection = getDesiredDirection(x1, x2, x3, x4, Lx, Ly, Lz);
    if dot(normv, desiredDirection)<0
        normv = -normv;
    end
    normv_all(elem,:) = normv;

    X = [x1(1), x2(1), x3(1), x4(1)];
    Y = [x1(2), x2(2), x3(2), x4(2)];
    Z = [x1(3), x2(3), x3(3), x4(3)];

    idx = 1;
    for ix=1:ngp_x
        epsilon = gp_x(ix);
        for iy=1:ngp_y
            eta = gp_y(iy);

            N1 = 0.25*(1 - epsilon)*(1 - eta);
            N2 = 0.25*(1 + epsilon)*(1 - eta);
            N3 = 0.25*(1 + epsilon)*(1 + eta);
            N4 = 0.25*(1 - epsilon)*(1 + eta);
            N = [N1, N2, N3, N4];

            x_g = N*X';
            y_g = N*Y';
            z_g = N*Z';

            x_gauss_all(elem, (idx-1)*3+1:(idx-1)*3+3) = [x_g, y_g, z_g];
            idx = idx+1;
        end
    end
end
end
