function velocity = calculateDomain(p, t, viscosity, x_col_inside, traction, u_known, normv_all, x_gauss_all, gp_x, gw_x, gp_y, gw_y, N, dN_dEpsilon, dN_dEta)

totalNodes = size(p,1);
totalElements = size(t,1);
numInnerPoints = size(x_col_inside,1);


AA = zeros(numInnerPoints*3, totalNodes*3);
BB = zeros(numInnerPoints*3, totalNodes*3);

[AA, BB] = buildDomainSystemMatrices(p, t, viscosity, x_col_inside, normv_all, x_gauss_all, gp_x, gw_x, gp_y, gw_y, N, dN_dEpsilon, dN_dEta, AA, BB, numInnerPoints, totalNodes, totalElements);

totalTraction = [traction(:,1); traction(:,2); traction(:,3)];
totalVelocity = [u_known(:,1); u_known(:,2); u_known(:,3)];
velocity = BB * totalTraction - AA * totalVelocity;

end

%% Helpers

function [AA, BB] = buildDomainSystemMatrices(p, t, viscosity, x_col_inside, normv_all, x_gauss_all, gp_x, gw_x, gp_y, gw_y, N, dN_dEpsilon, dN_dEta, AA, BB, numInnerPoints, totalNodes, totalElements)
    
    ngp_x = length(gp_x);
    ngp_y = length(gp_y);
    
    % Integral loop
    for domainPointIndex = 1:numInnerPoints
        x_col_r = x_col_inside(domainPointIndex,:);
        
        % Element loop
        for elementNum = 1:totalElements

            [normv, elemNodes, X, Y, Z] = getElementData(p, t, normv_all, elementNum);
            
            [T_ij_r, U_ij_r] = integrateKernelsForDomain(x_col_r, X, Y, Z, normv, viscosity, N, dN_dEpsilon, dN_dEta, gp_x, gw_x, gp_y, gw_y, ngp_x, ngp_y, x_gauss_all, elementNum);
            
            AA = assembleDomainMatrixA(AA, T_ij_r, elemNodes, domainPointIndex, totalNodes, numInnerPoints);
            BB = assembleDomainMatrixB(BB, U_ij_r, elemNodes, domainPointIndex, totalNodes, numInnerPoints);
        end
    end
end

function [normv, elemNodes, X, Y, Z] = getElementData(p, t, normv_all, elementNum)
    normv = normv_all(elementNum,:);
    elemNodes = t(elementNum,:);
    X = p(elemNodes,1);
    Y = p(elemNodes,2);
    Z = p(elemNodes,3);
end

function [T_ij_r, U_ij_r] = integrateKernelsForDomain(x_col_r, X, Y, Z, normv, viscosity, N, dN_dEpsilon, dN_dEta, gw_x, gw_y, ngp_x, ngp_y, x_gauss_all, elementNum)
    
    T_11_r=[0 0 0 0]; T_12_r=[0 0 0 0]; T_13_r=[0 0 0 0];
    T_21_r=[0 0 0 0]; T_22_r=[0 0 0 0]; T_23_r=[0 0 0 0];
    T_31_r=[0 0 0 0]; T_32_r=[0 0 0 0]; T_33_r=[0 0 0 0];

    U_11_r=[0 0 0 0]; U_12_r=[0 0 0 0]; U_13_r=[0 0 0 0];
    U_21_r=[0 0 0 0]; U_22_r=[0 0 0 0]; U_23_r=[0 0 0 0];
    U_31_r=[0 0 0 0]; U_32_r=[0 0 0 0]; U_33_r=[0 0 0 0];
    
    %% Gauss points loop
    for ix = 1:ngp_x
        w_x = gw_x(ix);
        for iy = 1:ngp_y
            w_y = gw_y(iy);
            idx = (ix-1)*ngp_y + iy;

            N_4 = N(idx,:);
            dN_dE = dN_dEpsilon(idx,:);
            dN_dH = dN_dEta(idx,:);
            [detJ, w] = calculateJacobianAndWeight(X, Y, Z, dN_dE, dN_dH, w_x, w_y);
            x_gauss_pt = getGaussPointCoordinates(x_gauss_all, elementNum, idx);
            
            Tij_r = T_ij3D(x_col_r, x_gauss_pt, normv);
            Uij_r = U_ij3D(x_col_r, x_gauss_pt, viscosity);
            
            [T_11_r, T_12_r, T_13_r, T_21_r, T_22_r, T_23_r, T_31_r, T_32_r, T_33_r, U_11_r, U_12_r, U_13_r, U_21_r, U_22_r, U_23_r, U_31_r, U_32_r, U_33_r] = calculateKernelContributions(T_11_r, T_12_r, T_13_r, T_21_r, T_22_r, T_23_r, T_31_r, T_32_r, T_33_r, U_11_r, U_12_r, U_13_r, U_21_r, U_22_r, U_23_r, U_31_r, U_32_r, U_33_r, Tij_r, Uij_r, w, detJ, N_4);
        end
    end
    
    T_ij_r = [T_11_r, T_12_r, T_13_r; T_21_r, T_22_r, T_23_r; T_31_r, T_32_r, T_33_r];
    U_ij_r = [U_11_r, U_12_r, U_13_r; U_21_r, U_22_r, U_23_r; U_31_r, U_32_r, U_33_r];
end

function [detJ, w] = calculateJacobianAndWeight(X, Y, Z, dN_dE, dN_dH, w_x, w_y)

    dx_dEpsilon = dN_dE*X; dy_dEpsilon = dN_dE*Y; dz_dEpsilon = dN_dE*Z;
    dx_dEta = dN_dH*X; dy_dEta = dN_dH*Y; dz_dEta = dN_dH*Z;
    
    T_epsilon = [dx_dEpsilon, dy_dEpsilon, dz_dEpsilon];
    T_eta = [dx_dEta, dy_dEta, dz_dEta];
    Nn = cross(T_epsilon, T_eta);
    detJ = norm(Nn);
    
    w = w_x * w_y;
end

function x_gauss_pt = getGaussPointCoordinates(x_gauss_all, elementNum, idx)
    x_gx = x_gauss_all(elementNum, (idx-1)*3+1);
    x_gy = x_gauss_all(elementNum, (idx-1)*3+2);
    x_gz = x_gauss_all(elementNum, (idx-1)*3+3);
    x_gauss_pt = [x_gx, x_gy, x_gz];
end

function [T_11_r, T_12_r, T_13_r, T_21_r, T_22_r, T_23_r, T_31_r, T_32_r, T_33_r, U_11_r, U_12_r, U_13_r, U_21_r, U_22_r, U_23_r, U_31_r, U_32_r, U_33_r] = calculateKernelContributions(T_11_r, T_12_r, T_13_r, T_21_r, T_22_r, T_23_r, T_31_r, T_32_r, T_33_r, U_11_r, U_12_r, U_13_r, U_21_r, U_22_r, U_23_r, U_31_r, U_32_r, U_33_r, Tij_r, Uij_r, w, detJ, N_4)
    for k = 1:4
        T_11_r(k) = T_11_r(k) + Tij_r(1,1)*w*detJ*N_4(k);
        T_12_r(k) = T_12_r(k) + Tij_r(1,2)*w*detJ*N_4(k);
        T_13_r(k) = T_13_r(k) + Tij_r(1,3)*w*detJ*N_4(k);
        T_21_r(k) = T_21_r(k) + Tij_r(2,1)*w*detJ*N_4(k);
        T_22_r(k) = T_22_r(k) + Tij_r(2,2)*w*detJ*N_4(k);
        T_23_r(k) = T_23_r(k) + Tij_r(2,3)*w*detJ*N_4(k);
        T_31_r(k) = T_31_r(k) + Tij_r(3,1)*w*detJ*N_4(k);
        T_32_r(k) = T_32_r(k) + Tij_r(3,2)*w*detJ*N_4(k);
        T_33_r(k) = T_33_r(k) + Tij_r(3,3)*w*detJ*N_4(k);
        

        U_11_r(k) = U_11_r(k) + Uij_r(1,1)*w*detJ*N_4(k);
        U_12_r(k) = U_12_r(k) + Uij_r(1,2)*w*detJ*N_4(k);
        U_13_r(k) = U_13_r(k) + Uij_r(1,3)*w*detJ*N_4(k);
        U_21_r(k) = U_21_r(k) + Uij_r(2,1)*w*detJ*N_4(k);
        U_22_r(k) = U_22_r(k) + Uij_r(2,2)*w*detJ*N_4(k);
        U_23_r(k) = U_23_r(k) + Uij_r(2,3)*w*detJ*N_4(k);
        U_31_r(k) = U_31_r(k) + Uij_r(3,1)*w*detJ*N_4(k);
        U_32_r(k) = U_32_r(k) + Uij_r(3,2)*w*detJ*N_4(k);
        U_33_r(k) = U_33_r(k) + Uij_r(3,3)*w*detJ*N_4(k);
    end
end

%% Will also be considered
function AA = assembleDomainMatrixA(AA, T_ij_r, elemNodes, domainPointIndex, totalNodes, numInnerPoints)
    
    AA(domainPointIndex,               elemNodes) = AA(domainPointIndex,               elemNodes) + T_ij_r(1,1);
    AA(domainPointIndex, totalNodes+elemNodes) = AA(domainPointIndex, totalNodes+elemNodes) +       T_ij_r(1,2);
    AA(domainPointIndex,2*totalNodes+elemNodes) = AA(domainPointIndex,2*totalNodes+elemNodes) +     T_ij_r(1,3);
    
    AA(domainPointIndex+numInnerPoints,          elemNodes) = AA(domainPointIndex+numInnerPoints,          elemNodes) +             T_ij_r(2,1);
    AA(domainPointIndex+numInnerPoints,totalNodes+elemNodes) = AA(domainPointIndex+numInnerPoints,totalNodes+elemNodes) +           T_ij_r(2,2);
    AA(domainPointIndex+numInnerPoints,2*totalNodes+elemNodes) = AA(domainPointIndex+numInnerPoints,2*totalNodes+elemNodes) +       T_ij_r(2,3);
    
    AA(domainPointIndex+2*numInnerPoints,          elemNodes) = AA(domainPointIndex+2*numInnerPoints,          elemNodes) +         T_ij_r(3,1);
    AA(domainPointIndex+2*numInnerPoints,totalNodes+elemNodes) = AA(domainPointIndex+2*numInnerPoints,totalNodes+elemNodes) +       T_ij_r(3,2);
    AA(domainPointIndex+2*numInnerPoints,2*totalNodes+elemNodes) = AA(domainPointIndex+2*numInnerPoints,2*totalNodes+elemNodes) +   T_ij_r(3,3);
end

function BB = assembleDomainMatrixB(BB, U_ij_r, elemNodes, domainPointIndex, totalNodes, numInnerPoints)
    
    BB(domainPointIndex,               elemNodes) = BB(domainPointIndex,               elemNodes) + U_ij_r(1,1);
    BB(domainPointIndex, totalNodes+elemNodes) = BB(domainPointIndex, totalNodes+elemNodes) + U_ij_r(1,2);
    BB(domainPointIndex,2*totalNodes+elemNodes) = BB(domainPointIndex,2*totalNodes+elemNodes) + U_ij_r(1,3);
    
    BB(domainPointIndex+numInnerPoints,          elemNodes) = BB(domainPointIndex+numInnerPoints,          elemNodes) + U_ij_r(2,1);
    BB(domainPointIndex+numInnerPoints,totalNodes+elemNodes) = BB(domainPointIndex+numInnerPoints,totalNodes+elemNodes) + U_ij_r(2,2);
    BB(domainPointIndex+numInnerPoints,2*totalNodes+elemNodes) = BB(domainPointIndex+numInnerPoints,2*totalNodes+elemNodes) + U_ij_r(2,3);
    
    BB(domainPointIndex+2*numInnerPoints,          elemNodes) = BB(domainPointIndex+2*numInnerPoints,          elemNodes) + U_ij_r(3,1);
    BB(domainPointIndex+2*numInnerPoints,totalNodes+elemNodes) = BB(domainPointIndex+2*numInnerPoints,totalNodes+elemNodes) + U_ij_r(3,2);
    BB(domainPointIndex+2*numInnerPoints,2*totalNodes+elemNodes) = BB(domainPointIndex+2*numInnerPoints,2*totalNodes+elemNodes) + U_ij_r(3,3);
end