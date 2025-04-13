clc; clear;

% Problem Parameters
Lx = 60; Ly = 30; Lz = 30;
Nx = 10; Ny = 5; Nz = 5;
dp_dx = -1;
viscosity = 1;
ngp_x = 5; ngp_y = 5;

bcType = 'dirichlet'; % Options: 'dirichlet', 'neumann', 'mixed'

%% Problem Specific Sizes
cubeDimensions = [5,5,5];
cubeCenter = [Lx/2, Ly/2, Lz/2];

%% Get Geometry (Looks Good)
[p,t] = getGeometry(Lx, Ly, Lz, Nx, Ny, Nz);

%% Get Points Inside (Looks Good)
x_col_inside = getPointsInside(Lx, Ly, Lz, 10,5,10,0.9);

%% Plot Break ☕
plotGeometryAndConnectivity(p, t);

%% Get BCs (Looks Good) (Check for more than 9 pts on inlet-outlet)
boundaryConditions = getBoundaryConditions(p, Nx, Ny, Nz, dp_dx, viscosity, Lx, Ly, Lz, bcType);

%% Get Gauss points and Shape Functions
[gp_x, gw_x, gp_y, gw_y]  = getGaussPoints(ngp_x, ngp_y);
[N,N_telles_BL,N_telles_BR,N_telles_TR,N_telles_TL,dN_dEpsilon, dN_dEta] = getShapeFunctions(gp_x, gp_y);

%% Precompute element normals and gauss points
tic
[normv_all, x_gauss_all] = precomputeElementData(p, t, gp_x, gp_y, Lx, Ly, Lz);
toc

%% Calculate boundaries
tic
[t_known, u_known, totalNodes, B, A] = calculateBoundaries(p, t, boundaryConditions, viscosity, normv_all, x_gauss_all, gp_x, gw_x, gp_y, gw_y, N,N_telles_BL, N_telles_BR, N_telles_TR, N_telles_TL, dN_dEpsilon, dN_dEta);
traction = t_known;
velocity = u_known;
toc

%% Calculate domain velocities
tic
domain = calculateDomain(p, t, viscosity, x_col_inside, traction, velocity, normv_all, x_gauss_all, gp_x, gw_x, gp_y, gw_y, N, dN_dEpsilon, dN_dEta);
toc

%% Plotting
tic
numInnerPoints = size(x_col_inside,1);
vx = domain(1:numInnerPoints);
vy = domain(numInnerPoints+1 : 2*numInnerPoints);
vz = domain(2*numInnerPoints+1 : 3*numInnerPoints);

figure;
hold on;
axis equal;
view(3);
xlabel('X'); ylabel('Y'); zlabel('Z');

quiver3(x_col_inside(:,1), x_col_inside(:,2), x_col_inside(:,3), vx, vy, vz, 'r');
grid on;
toc


a = Ly/2;
b = Lz/2;
cc = Lx/2;


for i = 1:Nx*Ny*Nz+2
    px = p(i,1);
    py = p(i,2);
    pz = p(i,3);

    % Shift coordinates to center
    p_x = px - cc;
    p_y = py - a;
    p_z = pz - b;
    uy = uy_analytical(p_y, p_z, a, b, dp_dx, viscosity);
    uz = uz_analytical(p_y, p_z, a, b, dp_dx, viscosity);


    if i <= Ny*Nz % Inlet face ( -p, visc uy, visc uz)
        fprintf("Inlet Node [%g, %g, %g]\n", px, py, px);
        t_calc = [0; viscosity*uy; viscosity*uz];
    elseif i <= Ny*Nz*2  % Outlet face ( p, -visc uy, -visc uz)
        fprintf("Outlet Node [%g, %g, %g]\n", px, py, px);
        t_calc = [0; -viscosity*uy; -viscosity*uz];
    elseif i <= Ny*Nz*2 + Ny * Nx % Bottom face ( -p, visc uy, visc uz)
        fprintf("Bottom Node [%g, %g, %g]\n", px, py, px);
        t_calc = [-viscosity*uy; 0; 0];
    elseif i <= Ny*Nz*2 + Ny*Nx*2
        % Top face
        fprintf("Top Node [%g, %g, %g]\n", px, py, px);
        t_calc = [viscosity*uy; 0; 0];
    elseif i <= Ny*Nz*2 + Ny*Nx*2 + Nx * Nz
        % Left face
        fprintf("Left Node [%g, %g, %g]\n", px, py, px);
        t_calc = [-viscosity*uz; 0; 0];
    else
        % Right face
        fprintf("Right Node [%g, %g, %g]\n", px, py, px);
        t_calc = [viscosity*uz; 0; 0];
    end

    t_numerical = traction([i, i+size(p,1), i+2*size(p,1)]);

    diff = t_numerical - t_calc;
    fprintf('Node %d:\n', i);
    fprintf('  Analytical t = [%g, %g, %g]\n', t_calc(1), t_calc(2), t_calc(3));
    fprintf('  Numerical t  = [%g, %g, %g]\n', t_numerical(1), t_numerical(2), t_numerical(3));
    fprintf('  Difference       = [%g, %g, %g]\n\n', diff(1), diff(2), diff(3));
end


function val = uy_analytical(p_y, p_z, a, b, dp_dx, viscosity)
    val = 0;
    factor = (16*a^2/(viscosity*pi^3))*(-dp_dx);
    for j = [1,3,5,7,9]
        sign_term = (-1)^((j-1)/2);
        term_z = 1 - cosh((j*pi*p_z/(2*a)))/cosh(j*pi*b/(2*a));
        dy_term = -sin(j*pi*p_y/(2*a))*(j*pi/(2*a));
        val = val + factor * sign_term * term_z * (dy_term / j^3);
    end
end

function val = uz_analytical(p_y, p_z, a, b, dp_dx, viscosity)

    val = 0;
    factor = (16*a^2/(viscosity*pi^3))*(-dp_dx);
    for j = [1,3,5,7,9]
        sign_term = (-1)^((j-1)/2);
        dterm_z_dz = -(sinh(j*pi*p_z/(2*a))*(j*pi/(2*a)))/cosh(j*pi*b/(2*a));
        term_y = cos(j*pi*p_y/(2*a));
        val = val + factor * sign_term * dterm_z_dz * term_y / j^3;
    end
end