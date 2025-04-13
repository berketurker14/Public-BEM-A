clc; clear;

%% Params
Lx = 60; Ly = 30; Lz = 30;
Nx = 10; Ny = 5; Nz = 5;
dp_dx = -1;
viscosity = 1;
ngp_x = 5; ngp_y = 5;
dimension = 3;
%% p and t
[p,t] = getGeometry(Lx, Ly, Lz, Nx, Ny, Nz);

%% Interior points
x_col_inside = getPointsInside(Lx, Ly, Lz, 10,1,10,0.9);

N = size(p,1);  % shorthand of Umut Hoca

%% Check plot
plotGeometryAndConnectivity(p,t);
%% {u_x, u_y, t_x, t_y} => 4*N unknowns

A_bem = zeros(dimension*N, dimension*2*N);  % (2N eqns, 4N unknowns)

%% BCs
boundaryConditions = getBoundaryConditions(p, Nx, Ny, Nz, dp_dx, viscosity, Lx, Ly, Lz);

%% Boundary System
tic;
for collNode = 1:N
    x_col = p(collNode,:);
    for elem = 1:N
        x1 = p(t(elem,1),:);
        x2 = p(t(elem,2),:);

        % Telles if same element
        if collNode == t(elem,1)
            x_gauss = N_telles_Bottom*[x1;x2];
            N_use   = N_telles_Bottom;
        elseif collNode == t(elem,2)
            x_gauss = N_telles_Top*[x1;x2];
            N_use   = N_telles_Top;
        else
            x_gauss = N_shape*[x1;x2];
            N_use   = N_shape;
        end

        pm   = x1 - x2;
        qlen = sqrt(pm(1)^2 + pm(2)^2);
        normv= [-pm(2), pm(1)]/qlen;

        T_11=[0 0]; T_12=[0 0]; T_21=[0 0]; T_22=[0 0];
        U_11=[0 0]; U_12=[0 0]; U_21=[0 0]; U_22=[0 0];

        % Integral Loop
        for j=1:ngp
            Tij= T_ij(x_col, x_gauss(j,:), normv);
            Uij= U_ij(x_col, x_gauss(j,:), viscosity);
            for k=1:2
                wght= gw(j)* J(elem)* N_use(j,k);
                T_11(1,k)= T_11(1,k)+ wght*Tij(1,1);
                T_12(1,k)= T_12(1,k)+ wght*Tij(1,2);
                T_21(1,k)= T_21(1,k)+ wght*Tij(2,1);
                T_22(1,k)= T_22(1,k)+ wght*Tij(2,2);

                U_11(1,k)= U_11(1,k)+ wght*Uij(1,1);
                U_12(1,k)= U_12(1,k)+ wght*Uij(1,2);
                U_21(1,k)= U_21(1,k)+ wght*Uij(2,1);
                U_22(1,k)= U_22(1,k)+ wght*Uij(2,2);
            end
        end

        rowUx = collNode;       
        rowUy = collNode + N;

        nd1 = t(elem,1);   nd2= t(elem,2);
        colUx1= nd1;       colUx2= nd2;
        colUy1= nd1 + N;   colUy2= nd2 + N;
        colTx1= nd1 + 2*N;  colTx2= nd2 + 2*N;
        colTy1= nd1 + 3*N;  colTy2= nd2 + 3*N;

        % T*u
        A_bem(rowUx, colUx1) = A_bem(rowUx, colUx1) + T_11(1,1);
        A_bem(rowUx, colUx2) = A_bem(rowUx, colUx2) + T_11(1,2);
        A_bem(rowUx, colUy1) = A_bem(rowUx, colUy1) + T_12(1,1);
        A_bem(rowUx, colUy2) = A_bem(rowUx, colUy2) + T_12(1,2);
        A_bem(rowUy, colUx1) = A_bem(rowUy, colUx1) + T_21(1,1);
        A_bem(rowUy, colUx2) = A_bem(rowUy, colUx2) + T_21(1,2);
        A_bem(rowUy, colUy1) = A_bem(rowUy, colUy1) + T_22(1,1);
        A_bem(rowUy, colUy2) = A_bem(rowUy, colUy2) + T_22(1,2);

        % -U*t => to be subtracted
        A_bem(rowUx, colTx1) = A_bem(rowUx, colTx1) - U_11(1,1);
        A_bem(rowUx, colTx2) = A_bem(rowUx, colTx2) - U_11(1,2);
        A_bem(rowUx, colTy1) = A_bem(rowUx, colTy1) - U_12(1,1);
        A_bem(rowUx, colTy2) = A_bem(rowUx, colTy2) - U_12(1,2);
        A_bem(rowUy, colTx1) = A_bem(rowUy, colTx1) - U_21(1,1);
        A_bem(rowUy, colTx2) = A_bem(rowUy, colTx2) - U_21(1,2);
        A_bem(rowUy, colTy1) = A_bem(rowUy, colTy1) - U_22(1,1);
        A_bem(rowUy, colTy2) = A_bem(rowUy, colTy2) - U_22(1,2);
    end
    A_bem(rowUx, collNode)     = A_bem(rowUx, collNode) + 0.5;    % + 1/2 C_ij part here
    A_bem(rowUy, collNode+N)   = A_bem(rowUy, collNode+N) + 0.5;  % + 1/2 C_ij part here
end
toc;
r_top = zeros(2*N,1);

%% E and c for lagrange
E = zeros(2*N, 4*N);
c = zeros(2*N,1);

for i=1:N
    yv = p(i,2);
    if i<=numOfNodes
        % bottom => no slip => (u_x=0,u_y=0)
        rowU = [zeros(1,4*N),0]; 
        % [1...2N -> ux,uy] [2N+1...4N -> tx,ty]
        r1 = zeros(1,4*N); r1(i)= 1;     % picks out u_x(i)
        r2 = zeros(1,4*N); r2(N + i)= 1; % picks out u_y(i)
        E=[E; r1; r2];
        c = [c; 0; 0];

    elseif i<=2*numOfNodes
        % right => parabolic => (u_x=4*(yv/h-(yv/h)^2),u_y=0)
        valx = 4*(yv/h - (yv/h)^2);
        r1 = zeros(1,4*N);  r1(i)=1;
        r2 = zeros(1,4*N);  r2(N + i)=1;
        E=[E; r1; r2];
        c  = [c; valx; 0];

    elseif i<=3*numOfNodes
        % top => (u_x=0,u_y=0)
        r1=zeros(1,4*N); r1(i)=1; 
        r2=zeros(1,4*N); r2(N+i)=1;
        E=[E; r1; r2];
        c=[c; 0; 0];

    elseif i<=4*numOfNodes
        % left => parabolic => (u_x=..., u_y=0)
        valx = 4*(yv/h-(yv/h)^2);
        r1=zeros(1,4*N); r1(i)=1; 
        r2=zeros(1,4*N); r2(N+i)=1;
        E=[E; r1; r2];
        c=[c; valx; 0];

    else
        % cylinder => no slip => (u=0)
        r1=zeros(1,4*N); r1(i)=1;
        r2=zeros(1,4*N); r2(N+i)=1;
        E=[E; r1; r2];
        c=[c; 0; 0];
    end
end



%% Lagrange multipliers
numConstraints = size(E,1);
Z_topright = zeros(2*N, numConstraints);    % 2N x 2N
Z_bottomright = zeros(numConstraints, numConstraints);

bigMatrix = [
    A_bem,     Z_topright;      % 2N x (4N + 2N) = 2N x 6N
    E,         Z_bottomright    % 2N x (4N + 2N) = 2N x 6N
];
bigRHS = [ zeros(2*N,1); c];


%% Solve
solution = bigMatrix \ bigRHS;

u_and_t = solution(1:4*N);
lambda  = solution(4*N+1:end);

u_part = u_and_t(1:2*N);
t_part = u_and_t(2*N+1:4*N);

%% Interior System
AA = zeros(2*numInnerPoints, 2*N);
BB = zeros(2*numInnerPoints, 2*N);

for domainPt=1:numInnerPoints
    x_col_r = x_col_inside(domainPt,:);
    for elem=1:N
        x1 = p(t(elem,1),:);
        x2 = p(t(elem,2),:);
        pm = x1 - x2;
        ql = sqrt(pm(1)^2+pm(2)^2);
        normv= [-pm(2), pm(1)]/ql;

        T_11=[0 0]; T_12=[0 0]; T_21=[0 0]; T_22=[0 0];
        U_11=[0 0]; U_12=[0 0]; U_21=[0 0]; U_22=[0 0];

        for j=1:ngp
            Xg = N_shape(j,:)*[x1;x2];
            Tij= T_ij(x_col_r, Xg, normv);
            Uij= U_ij(x_col_r, Xg, viscosity);
            for k=1:2
                wght= gw(j)*J(elem)*N_shape(j,k);
                T_11(1,k)= T_11(1,k)+ wght*Tij(1,1);
                T_12(1,k)= T_12(1,k)+ wght*Tij(1,2);
                T_21(1,k)= T_21(1,k)+ wght*Tij(2,1);
                T_22(1,k)= T_22(1,k)+ wght*Tij(2,2);

                U_11(1,k)= U_11(1,k)+ wght*Uij(1,1);
                U_12(1,k)= U_12(1,k)+ wght*Uij(1,2);
                U_21(1,k)= U_21(1,k)+ wght*Uij(2,1);
                U_22(1,k)= U_22(1,k)+ wght*Uij(2,2);
            end
        end

        rowUx = domainPt;
        rowUy = domainPt + numInnerPoints;

        nd1 = t(elem,1);
        nd2 = t(elem,2);

        colUx1=nd1; colUx2=nd2;
        colUy1=nd1+totalNodes; colUy2=nd2+totalNodes;

        AA(rowUx,colUx1)= AA(rowUx,colUx1)+ T_11(1,1);
        AA(rowUx,colUx2)= AA(rowUx,colUx2)+ T_11(1,2);
        AA(rowUx,colUy1)= AA(rowUx,colUy1)+ T_12(1,1);
        AA(rowUx,colUy2)= AA(rowUx,colUy2)+ T_12(1,2);
        AA(rowUy,colUx1)= AA(rowUy,colUx1)+ T_21(1,1);
        AA(rowUy,colUx2)= AA(rowUy,colUx2)+ T_21(1,2);
        AA(rowUy,colUy1)= AA(rowUy,colUy1)+ T_22(1,1);
        AA(rowUy,colUy2)= AA(rowUy,colUy2)+ T_22(1,2);

        BB(rowUx,colUx1)= BB(rowUx,colUx1)+ U_11(1,1);
        BB(rowUx,colUx2)= BB(rowUx,colUx2)+ U_11(1,2);
        BB(rowUx,colUy1)= BB(rowUx,colUy1)+ U_12(1,1);
        BB(rowUx,colUy2)= BB(rowUx,colUy2)+ U_12(1,2);
        BB(rowUy,colUx1)= BB(rowUy,colUx1)+ U_21(1,1);
        BB(rowUy,colUx2)= BB(rowUy,colUx2)+ U_21(1,2);
        BB(rowUy,colUy1)= BB(rowUy,colUy1)+ U_22(1,1);
        BB(rowUy,colUy2)= BB(rowUy,colUy2)+ U_22(1,2);
    end
end

velocity = BB*t_part - AA*u_part;

u_x = velocity(1:numInnerPoints);
u_y = velocity(numInnerPoints+1:2*numInnerPoints);

%% Plotting
figure;
quiver(x_col_inside(:,1), x_col_inside(:,2), u_x, u_y, 'r');
hold on; 
plot(p(:,1), p(:,2), 'bo');
axis equal; grid on;
title('Velocity Field');
xlabel('X'); ylabel('Y');