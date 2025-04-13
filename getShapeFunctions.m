function [N,N_telles_BL,N_telles_BR,N_telles_TR,N_telles_TL,dN_dEpsilon, dN_dEta] = getShapeFunctions(gp_x, gp_y)

ngp_x = length(gp_x);
ngp_y = length(gp_y);

N = zeros(ngp_x*ngp_y,4);
N_telles_BL = zeros(ngp_x*ngp_y,4);
N_telles_BR = zeros(ngp_x*ngp_y,4);
N_telles_TR = zeros(ngp_x*ngp_y,4);
N_telles_TL = zeros(ngp_x*ngp_y,4);
% Derivatives are required for jacobian
dN_dEpsilon = zeros(ngp_x*ngp_y,4);
dN_dEta = zeros(ngp_x*ngp_y,4);

idx=1;

telles_gp_y_minus1 = tellesTransformation(-1,gp_y);
telles_gp_x_minus1 = tellesTransformation(-1,gp_x);
telles_gp_y_plus1 = tellesTransformation(1,gp_y);
telles_gp_x_plus1 = tellesTransformation(1,gp_x);
for ix=1:ngp_x
    epsilon=gp_x(ix);
    telles_epsilon_BL = telles_gp_x_minus1(ix);
    telles_epsilon_BR = telles_gp_x_plus1(ix);
    telles_epsilon_TR = telles_gp_x_plus1(ix);
    telles_epsilon_TL = telles_gp_x_minus1(ix);
    for iy=1:ngp_y
        telles_eta_BL = telles_gp_y_minus1(iy);
        telles_eta_BR = telles_gp_y_minus1(iy);
        telles_eta_TR = telles_gp_y_plus1(iy);
        telles_eta_TL = telles_gp_y_plus1(iy);
        eta=gp_y(iy);
        N1 = 0.25*(1 - epsilon)*(1 - eta);
        N2 = 0.25*(1 + epsilon)*(1 - eta);
        N3 = 0.25*(1 + epsilon)*(1 + eta);
        N4 = 0.25*(1 - epsilon)*(1 + eta);
        N_telles_BL(idx,:) = [0.25*(1 - telles_epsilon_BL)*(1 - telles_eta_BL);0.25*(1 + telles_epsilon_BL)*(1 - telles_eta_BL);0.25*(1 + telles_epsilon_BL)*(1 + telles_eta_BL);0.25*(1 - telles_epsilon_BL)*(1 + telles_eta_BL);];
        N_telles_BR(idx,:) = [0.25*(1 - telles_epsilon_BR)*(1 - telles_eta_BR);0.25*(1 + telles_epsilon_BR)*(1 - telles_eta_BR);0.25*(1 + telles_epsilon_BR)*(1 + telles_eta_BR);0.25*(1 - telles_epsilon_BR)*(1 + telles_eta_BR);];
        N_telles_TR(idx,:) = [0.25*(1 - telles_epsilon_TR)*(1 - telles_eta_TR);0.25*(1 + telles_epsilon_TR)*(1 - telles_eta_TR);0.25*(1 + telles_epsilon_TR)*(1 + telles_eta_TR);0.25*(1 - telles_epsilon_TR)*(1 + telles_eta_TR);];
        N_telles_TL(idx,:) = [0.25*(1 - telles_epsilon_TL)*(1 - telles_eta_TL);0.25*(1 + telles_epsilon_TL)*(1 - telles_eta_TL);0.25*(1 + telles_epsilon_TL)*(1 + telles_eta_TL);0.25*(1 - telles_epsilon_TL)*(1 + telles_eta_TL);];
        N(idx,:) = [N1, N2, N3, N4];

        dN_dEpsilon(idx,:) = 0.25*[-(1 - eta),(1 - eta),(1 + eta),-(1 + eta)];
        dN_dEta(idx,:) = 0.25*[-(1 - epsilon),-(1 + epsilon),(1 + epsilon),(1 - epsilon)];

        idx=idx+1;
    end
end
end
