function [t_known, u_known, totalNodes, B, A] = calculateBoundaries(p, t, boundaryConditions, viscosity, normv_all, x_gauss_all, gp_x, gw_x, gp_y, gw_y, N, N_telles_BL, N_telles_BR, N_telles_TR, N_telles_TL, dN_dEpsilon, dN_dEta)

totalNodes = size(p,1);
totalElements = size(t,1);

A = zeros(totalNodes * 3);
B = zeros(totalNodes * 3);
u_known = boundaryConditions(:,1:3);
t_known = boundaryConditions(:,4:6);

% Extract known values and their indices
[u_known_idx, t_known_idx, u_count, t_count] = extractBoundaryIndices(u_known, t_known, totalNodes);

% Build the system matrices
[A, B] = buildSystemMatrices(p, t, u_known, t_known, u_known_idx, t_known_idx, u_count, t_count,viscosity, normv_all, N, N_telles_BL, N_telles_BR, N_telles_TR, N_telles_TL, dN_dEpsilon, dN_dEta, gp_x, gw_x, gp_y, gw_y, totalNodes, totalElements);

% Solve the system for unknown values
[u_known, t_known] = solveSystemForUnknowns(A, B, u_known, t_known, u_known_idx, t_known_idx, u_count, t_count, totalNodes);

% Return velocity and traction as separate arrays
velocity = u_known;
traction = t_known;

end

%% Helper Functions

function [u_known_idx, t_known_idx, u_count, t_count] = extractBoundaryIndices(u_known, t_known, totalNodes)
    % Extract indices of known values from the boundary conditions
    
    t_count = 0;
    u_count = 0;
    u_known_idx = cell(totalNodes, 1);
    t_known_idx = cell(totalNodes, 1);
    
    % Find indices where velocities and tractions are known
    for i = 1:totalNodes
        cols_not_nan = find(~isnan(u_known(i,:)));    
        if ~isempty(cols_not_nan)
            u_count = u_count+1;
            u_known_idx{u_count,1} = [cols_not_nan,i];
        end
        
        cols_not_nan = find(~isnan(t_known(i,:)));    
        if ~isempty(cols_not_nan)
            t_count = t_count+1;        
            t_known_idx{t_count,1} = [cols_not_nan,i];
        end
    end
end

function [A, B] = buildSystemMatrices(p, t, u_known, t_known, u_known_idx, t_known_idx, u_count, t_count, viscosity, normv_all, N, N_telles_BL, N_telles_BR, N_telles_TR, N_telles_TL, dN_dEpsilon, dN_dEta, gp_x, gw_x, gp_y, gw_y, totalNodes, totalElements)
    % Build the BEM system matrices A(multiplier of u) and B(multiplier of t)
    
    A = zeros(totalNodes * 3);
    B = zeros(totalNodes * 3);
    ngp_x = length(gp_x);
    ngp_y = length(gp_y);

    % Precompute bc info
    bc_info = precomputeBoundaryInfo(u_known_idx, t_known_idx, u_known, t_known, u_count, t_count);
    
    %% Collocation node loop
    for collocation_node = 1:totalNodes
        x_col = p(collocation_node, :);
        
        %% Element loop
        for elementNum = 1:totalElements
            normv = normv_all(elementNum,:);
            elemNodes = t(elementNum,:);
            X = p(elemNodes,1); 
            Y = p(elemNodes,2);
            Z = p(elemNodes,3);
            
            % Get bc info to use in element.
            [ux, uy, uz, tx, ty, tz] = getBoundaryValuesForElement(elemNodes, bc_info, u_known_idx, t_known_idx);
            
            % Integrate over the element
            [T_ij, U_ij] = integrateOverElement(collocation_node, elementNum, x_col, X, Y, Z, normv, viscosity, N, N_telles_BL, N_telles_BR, N_telles_TR, N_telles_TL, dN_dEpsilon, dN_dEta, gp_x, gw_x, gp_y, gw_y, ngp_x, ngp_y, ux, uy, uz, tx, ty, tz, t);
            
            % Assemble into global system matrices
            A = assembleMatrixA(A, T_ij, elemNodes, collocation_node, totalNodes);
            B = assembleMatrixB(B, U_ij, elemNodes, collocation_node, totalNodes);
        end
    end
end

function bc_info = precomputeBoundaryInfo(u_known_idx, t_known_idx, u_known, t_known, u_count, t_count)
    
    bc_info = struct();
    u_known_node_indices = zeros(u_count, 1);
    t_known_node_indices = zeros(t_count, 1);
    
    for i = 1:u_count
        u_known_node_indices(i) = u_known_idx{i}(end);
        bc_info.u_components{i} = u_known_idx{i}(1:end-1); % Components that are known
        bc_info.u_values{i} = u_known(u_known_idx{i}(end), :); % Velocity values
    end
    
    for i = 1:t_count
        t_known_node_indices(i) = t_known_idx{i}(end);
        bc_info.t_components{i} = t_known_idx{i}(1:end-1); % Components that are known
        bc_info.t_values{i} = t_known(t_known_idx{i}(end), :); % Traction values
    end
    
    bc_info.u_known_node_indices = u_known_node_indices;
    bc_info.t_known_node_indices = t_known_node_indices;
end

function [ux, uy, uz, tx, ty, tz] = getBoundaryValuesForElement(elemNodes, bc_info, u_known_idx, t_known_idx)

    % Default values - multipliers for kernels
    ux = 1; uy = 1; uz = 1; tx = 1; ty = 1; tz = 1;
    
    % First check for velocity boundary conditions
    for i = 1:length(elemNodes)
        nodeId = elemNodes(i);
        % Find if this node has known velocity components
        node_idx = find(bc_info.u_known_node_indices == nodeId);
        if ~isempty(node_idx)
            % Set known velocity components
            comps = bc_info.u_components{node_idx};
            values = bc_info.u_values{node_idx};
            if ismember(1, comps)
                ux = values(1);
            end
            if ismember(2, comps)
                uy = values(2);
            end
            if ismember(3, comps)
                uz = values(3);
            end
        end
    end
    
    % Then check for traction boundary conditions
    for i = 1:length(elemNodes)
        nodeId = elemNodes(i);
        % Find if this node has known traction components
        node_idx = find(bc_info.t_known_node_indices == nodeId);
        if ~isempty(node_idx)
            % Set known traction components
            comps = bc_info.t_components{node_idx};
            values = bc_info.t_values{node_idx};
            if ismember(1, comps)
                tx = values(1);
            end
            if ismember(2, comps)
                ty = values(2);
            end
            if ismember(3, comps)
                tz = values(3);
            end
        end
    end
end

function [T_ij, U_ij] = integrateOverElement(collocation_node, elementNum, x_col, X, Y, Z, normv, viscosity, N, N_telles_BL, N_telles_BR, N_telles_TR, N_telles_TL, dN_dEpsilon, dN_dEta, gp_x, gw_x, gp_y, gw_y, ngp_x, ngp_y, ux, uy, uz, tx, ty, tz, t)
    % Integrate kernel functions over an element
    
    % Initialize matrices for T and U kernels
    T_11=[0 0 0 0]; T_12=[0 0 0 0]; T_13=[0 0 0 0];
    T_21=[0 0 0 0]; T_22=[0 0 0 0]; T_23=[0 0 0 0];
    T_31=[0 0 0 0]; T_32=[0 0 0 0]; T_33=[0 0 0 0];

    U_11=[0 0 0 0]; U_12=[0 0 0 0]; U_13=[0 0 0 0];
    U_21=[0 0 0 0]; U_22=[0 0 0 0]; U_23=[0 0 0 0];
    U_31=[0 0 0 0]; U_32=[0 0 0 0]; U_33=[0 0 0 0];
    
    % Loop over Gauss points
    for ix = 1:ngp_x
        w_x = gw_x(ix);
        for iy = 1:ngp_y
            w_y = gw_y(iy);
            idx = (ix-1)*ngp_y + iy;
            
            % Apply Telles transformation for singular elements
            if collocation_node == t(elementNum,1) % bottom left is singular
                x_gauss_pt = N_telles_BL(idx,:) * [X, Y, Z];
                N_4 = N_telles_BL(idx,:);
            elseif collocation_node == t(elementNum,2) % bottom right is singular
                x_gauss_pt = N_telles_BR(idx,:) * [X, Y, Z];
                N_4 = N_telles_BR(idx,:);
            elseif collocation_node == t(elementNum,3) % top right is singular
                x_gauss_pt = N_telles_TR(idx,:) * [X, Y, Z];
                N_4 = N_telles_TR(idx,:);
            elseif collocation_node == t(elementNum,4) % top left is singular
                x_gauss_pt = N_telles_TL(idx,:) * [X, Y, Z];
                N_4 = N_telles_TL(idx,:);
            else
                x_gauss_pt = N(idx,:) * [X, Y, Z];
                N_4 = N(idx,:);
            end
            
            % Calculate Jacobian
            dN_dE = dN_dEpsilon(idx,:);
            dN_dH = dN_dEta(idx,:);
            
            dx_dEpsilon = dN_dE*X;
            dy_dEpsilon = dN_dE*Y;
            dz_dEpsilon = dN_dE*Z;

            dx_dEta = dN_dH*X;
            dy_dEta = dN_dH*Y;
            dz_dEta = dN_dH*Z;

            T_epsilon = [dx_dEpsilon, dy_dEpsilon, dz_dEpsilon];
            T_eta = [dx_dEta, dy_dEta, dz_dEta];
            Nn = cross(T_epsilon, T_eta);
            detJ = norm(Nn);
            
            % Weight for numerical integration
            w = w_x * w_y;
            
            % Evaluate kernel functions
            Tij = T_ij3D(x_col, x_gauss_pt, normv);
            Uij = U_ij3D(x_col, x_gauss_pt, viscosity);
            
            % Accumulate contributions for each shape function
            for k = 1:4
                % T kernel contributions
                T_11(k) = T_11(k) + Tij(1,1)*w*detJ*N_4(k)*ux;
                T_12(k) = T_12(k) + Tij(1,2)*w*detJ*N_4(k)*uy;
                T_13(k) = T_13(k) + Tij(1,3)*w*detJ*N_4(k)*uz;
                T_21(k) = T_21(k) + Tij(2,1)*w*detJ*N_4(k)*ux;
                T_22(k) = T_22(k) + Tij(2,2)*w*detJ*N_4(k)*uy;
                T_23(k) = T_23(k) + Tij(2,3)*w*detJ*N_4(k)*uz;
                T_31(k) = T_31(k) + Tij(3,1)*w*detJ*N_4(k)*ux;
                T_32(k) = T_32(k) + Tij(3,2)*w*detJ*N_4(k)*uy;
                T_33(k) = T_33(k) + Tij(3,3)*w*detJ*N_4(k)*uz;

                % U kernel contributions
                U_11(k) = U_11(k) + Uij(1,1)*w*detJ*N_4(k)*tx;
                U_12(k) = U_12(k) + Uij(1,2)*w*detJ*N_4(k)*ty;
                U_13(k) = U_13(k) + Uij(1,3)*w*detJ*N_4(k)*tz;
                U_21(k) = U_21(k) + Uij(2,1)*w*detJ*N_4(k)*tx;
                U_22(k) = U_22(k) + Uij(2,2)*w*detJ*N_4(k)*ty;
                U_23(k) = U_23(k) + Uij(2,3)*w*detJ*N_4(k)*tz;
                U_31(k) = U_31(k) + Uij(3,1)*w*detJ*N_4(k)*tx;
                U_32(k) = U_32(k) + Uij(3,2)*w*detJ*N_4(k)*ty;
                U_33(k) = U_33(k) + Uij(3,3)*w*detJ*N_4(k)*tz;
            end
        end
    end
    
    % Pack results into structure
    T_ij = [T_11, T_12, T_13; T_21, T_22, T_23; T_31, T_32,T_33];
              
    U_ij = [U_11, U_12, U_13; U_21, U_22, U_23; U_31, U_32,U_33];
end

%% C++ --> assembleMatrix on manager for assembleMatrixA and assembleMatrixB
function A = assembleMatrixA(A, T_ij, elemNodes, collocation_node, totalNodes)
    % T kernel contributions into A
    
    A(collocation_node,                           elemNodes) = A(collocation_node,                           elemNodes) +   T_ij(1,1);
    A(collocation_node,               totalNodes + elemNodes) = A(collocation_node,               totalNodes + elemNodes) + T_ij(1,2);
    A(collocation_node,             2*totalNodes + elemNodes) = A(collocation_node,             2*totalNodes + elemNodes) + T_ij(1,3);
    
    A(collocation_node+totalNodes,                elemNodes) = A(collocation_node+totalNodes,                elemNodes)   + T_ij(2,1);
    A(collocation_node+totalNodes,    totalNodes + elemNodes) = A(collocation_node+totalNodes,    totalNodes + elemNodes) + T_ij(2,2);
    A(collocation_node+totalNodes,  2*totalNodes + elemNodes) = A(collocation_node+totalNodes,  2*totalNodes + elemNodes) + T_ij(2,3);
    
    A(collocation_node+2*totalNodes,              elemNodes) = A(collocation_node+2*totalNodes,              elemNodes)   + T_ij(3,1);
    A(collocation_node+2*totalNodes,  totalNodes + elemNodes) = A(collocation_node+2*totalNodes,  totalNodes + elemNodes) + T_ij(3,2);
    A(collocation_node+2*totalNodes,2*totalNodes + elemNodes) = A(collocation_node+2*totalNodes,2*totalNodes + elemNodes) + T_ij(3,3);
end

function B = assembleMatrixB(B, U_ij, elemNodes, collocation_node, totalNodes)
    % U kernel contributions into B
    
    B(collocation_node,                           elemNodes) = B(collocation_node,                           elemNodes)   + U_ij(1,1);
    B(collocation_node,              totalNodes  + elemNodes) = B(collocation_node,              totalNodes  + elemNodes) + U_ij(1,2);
    B(collocation_node,             2*totalNodes + elemNodes) = B(collocation_node,             2*totalNodes + elemNodes) + U_ij(1,3);
    
    B(collocation_node+totalNodes,                elemNodes) = B(collocation_node+totalNodes,                elemNodes)   + U_ij(2,1);
    B(collocation_node+totalNodes,    totalNodes + elemNodes) = B(collocation_node+totalNodes,    totalNodes + elemNodes) + U_ij(2,2);
    B(collocation_node+totalNodes,  2*totalNodes + elemNodes) = B(collocation_node+totalNodes,  2*totalNodes + elemNodes) + U_ij(2,3);
    
    B(collocation_node+2*totalNodes,              elemNodes) = B(collocation_node+2*totalNodes,              elemNodes)   + U_ij(3,1);
    B(collocation_node+2*totalNodes,  totalNodes + elemNodes) = B(collocation_node+2*totalNodes,  totalNodes + elemNodes) + U_ij(3,2);
    B(collocation_node+2*totalNodes,2*totalNodes + elemNodes) = B(collocation_node+2*totalNodes,2*totalNodes + elemNodes) + U_ij(3,3);
end


function [u_known, t_known] = solveSystemForUnknowns(A, B, u_known, t_known, u_known_idx, t_known_idx, u_count, t_count, totalNodes)
    
    u_x_known = [];
    u_y_known = [];
    u_z_known = [];
    t_x_known = [];
    t_y_known = [];
    t_z_known = [];
    
    % Map indices to whole system
    % returns
    for i = 1:u_count
        node_idx = u_known_idx{i}(end);
        components = u_known_idx{i}(1:end-1);
        
        for j = 1:length(components)
            comp = components(j);
            if comp == 1
                u_x_known = [u_x_known; node_idx];
            elseif comp == 2
                u_y_known = [u_y_known; node_idx];
            elseif comp == 3
                u_z_known = [u_z_known; node_idx];
            end
        end
    end
    
    for i = 1:t_count
        node_idx = t_known_idx{i}(end);
        components = t_known_idx{i}(1:end-1);
        
        for j = 1:length(components)
            comp = components(j);
            if comp == 1
                t_x_known = [t_x_known; node_idx];
            elseif comp == 2
                t_y_known = [t_y_known; node_idx];
            elseif comp == 3
                t_z_known = [t_z_known; node_idx];
            end
        end
    end
    
    %% X-component
    u_indices = u_x_known;
    t_indices = t_x_known;
    
    % no initial bc u_x of node given
    if ~isempty(u_indices)
        A_ux = A(u_indices, u_indices);
        B_ux = B(u_indices, u_indices);
        u_vals = u_known(u_indices, 1);
        
        C = 0.5 * eye(length(u_indices));
        
        % Solve for x-traction
        rhs = (A_ux + C) * u_vals;
        t_vals =  B_ux \ rhs; % will be changed
        
        % Update traction values
        for i = 1:length(u_indices)
            node_idx = u_indices(i);
            t_known(node_idx, 1) = t_vals(i);
        end
    end
    
    if ~isempty(t_indices)
        % Extract system for x-traction known
        A_tx = A(t_indices, t_indices);
        B_tx = B(t_indices, t_indices);
        t_vals = t_known(t_indices, 1);
        
        % Identity matrix for jump term
        C = 0.5 * eye(length(t_indices));
        
        % Solve for x-velocity
        rhs = B_tx * t_vals;
        u_vals = (A_tx + C) \ rhs;
        
        % Update velocity values
        for i = 1:length(t_indices)
            node_idx = t_indices(i);
            u_known(node_idx, 1) = u_vals(i);
        end
    end
    
    %% Y-component
    u_indices = u_y_known;
    t_indices = t_y_known;
    
    if ~isempty(u_indices)
        % Extract system for y-velocity known
        A_uy = A(totalNodes + u_indices, totalNodes + u_indices);
        B_uy = B(totalNodes + u_indices, totalNodes + u_indices);
        u_vals = u_known(u_indices, 2);
        
        % Identity matrix for jump term
        C = 0.5 * eye(length(u_indices));
        
        % Solve for y-traction
        rhs = (A_uy + C) * u_vals;
        t_vals = B_uy \ rhs;
        
        % Update traction values
        for i = 1:length(u_indices)
            node_idx = u_indices(i);
            t_known(node_idx, 2) = t_vals(i);
        end
    end
    
    if ~isempty(t_indices)
        % Extract system for y-traction known
        A_ty = A(totalNodes + t_indices, totalNodes + t_indices);
        B_ty = B(totalNodes + t_indices, totalNodes + t_indices);
        t_vals = t_known(t_indices, 2);
        
        % Identity matrix for jump term
        C = 0.5 * eye(length(t_indices));
        
        % Solve for y-velocity
        rhs = B_ty * t_vals;
        u_vals = (A_ty + C) \ rhs;
        
        % Update velocity values
        for i = 1:length(t_indices)
            node_idx = t_indices(i);
            u_known(node_idx, 2) = u_vals(i);
        end
    end
    
    %% Z-component
    u_indices = u_z_known;
    t_indices = t_z_known;
    
    if ~isempty(u_indices)
        % Extract system for z-velocity known
        A_uz = A(2*totalNodes + u_indices, 2*totalNodes + u_indices);
        B_uz = B(2*totalNodes + u_indices, 2*totalNodes + u_indices);
        u_vals = u_known(u_indices, 3);
        
        % Identity matrix for jump term
        C = 0.5 * eye(length(u_indices));
        
        % Solve for z-traction
        rhs = (A_uz + C) * u_vals;
        t_vals = B_uz \ rhs;
        
        % Update traction values
        for i = 1:length(u_indices)
            node_idx = u_indices(i);
            t_known(node_idx, 3) = t_vals(i);
        end
    end
    
    if ~isempty(t_indices)
        % Extract system for z-traction known
        A_tz = A(2*totalNodes + t_indices, 2*totalNodes + t_indices);
        B_tz = B(2*totalNodes + t_indices, 2*totalNodes + t_indices);
        t_vals = t_known(t_indices, 3);
        
        % Identity matrix for jump term
        C = 0.5 * eye(length(t_indices));
        
        % Solve for z-velocity
        rhs = B_tz * t_vals;
        u_vals = (A_tz + C) \ rhs;
        
        % Update velocity values
        for i = 1:length(t_indices)
            node_idx = t_indices(i);
            u_known(node_idx, 3) = u_vals(i);
        end
    end
end