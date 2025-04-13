function boundaryConditions = getBoundaryConditions(p, Nx, Ny, Nz, dp_dx, viscosity, Lx, Ly, Lz, bcType)
    % bcType: 'dirichlet', 'neumann', or 'mixed'
    
    totalNodes = size(p,1);
    boundaryConditions = zeros(totalNodes, 6);
    boundaryConditions(:,:) = NaN;
    
    % Geometrical constants for frank m white square duct
    a = Ly/2;
    b = Lz/2;
    
    % Identify boundary regions
    inletNodes = 1 : (Nz+1)*(Ny+1);
    outletNodes = (Nz+1)*(Ny+1)+1 : (Nz+1)*(Ny+1)*2;
    rightNodes = (Nz+1)*(Ny+1)*2+1 : (Nz+1)*(Nx-1)+(Nz+1)*(Ny+1)*2;
    leftNodes = (Nz+1)*(Nx-1)+(Nz+1)*(Ny+1)*2+1 : 2*(Nz+1)*(Nx-1)+(Nz+1)*(Ny+1)*2;
    bottomNodes = 2*(Nz+1)*(Nx-1)+(Nz+1)*(Ny+1)*2+1 : 2*(Nz+1)*(Nx-1)+(Nz+1)*(Ny+1)*2+(Nx-1)*(Ny-1); % Split to left and right!!!!!!!!!!!!!!!
    topNodes = 2*(Nz+1)*(Nx-1)+(Nz+1)*(Ny+1)*2+(Nx-1)*(Ny-1)+1 : totalNodes;

    
    %% Apply boundary conditions based on type
    if strcmpi(bcType, 'dirichlet')
        % All Dirichlet - all velocities known
        boundaryConditions = applyDirichletBoundaries(boundaryConditions, p, a, b, dp_dx, viscosity, inletNodes, outletNodes);
    elseif strcmpi(bcType, 'neumann')
        % All Neumann - all tractions known
        boundaryConditions = applyNeumannBoundaries(boundaryConditions, p, a, b, dp_dx, viscosity, inletNodes, outletNodes, rightNodes, leftNodes, bottomNodes, topNodes);
    else
        % Mixed BCs
        boundaryConditions = applyMixedBoundaries(boundaryConditions, p, a, b, dp_dx, viscosity, inletNodes, outletNodes, rightNodes, leftNodes, bottomNodes, topNodes);
    end
end

function boundaryConditions = applyDirichletBoundaries(boundaryConditions, p, a, b, dp_dx, viscosity, inletNodes, outletNodes)
    % Set velocity (Dirichlet) for all boundaries Frank M. White
    constant = 16*a^2/(viscosity*pi^3)*(-dp_dx);
    
    % Apply to all nodes
    for i = 1:size(p,1)
        boundaryConditions(i,2) = 0; % v=0
        boundaryConditions(i,3) = 0; % w=0
        
        p_y = p(i,2) - a;
        p_z = p(i,3) - b;
        
        % Determine u velocity
        if ismember(i, [inletNodes, outletNodes])
            % Inlet/Outlet: use analytical profile
            total = 0;
            for j = 1:2:9
                total = total + constant * (-1)^((j-1)/2) * (1-(cosh(j*pi*p_z/(2*a))/cosh(j*pi*b/(2*a)))) * cos(j*pi*p_y/(2*a))/j^3;
            end
            boundaryConditions(i,1) = total;
        else
            % Other walls: no-slip
            boundaryConditions(i,1) = 0;
        end
    end
end

% Not functional yet, requires pressure to be known additional to velocity
% profile or tractions
function boundaryConditions = applyNeumannBoundaries(boundaryConditions, p, a, b, dp_dx, viscosity, inletNodes, outletNodes, bottomNodes, topNodes, leftNodes, rightNodes)
    % Set traction for all boundaries
    constant = 16*a^2/(viscosity*pi^3)*(-dp_dx);
    
    for i = 1:size(p,1)
        p_y = p(i,2) - a;
        p_z = p(i,3) - b;
        
        % Set all normal and tangential tractions based on location
        if ismember(i, inletNodes)
            % Inlet: pressure and shear
            boundaryConditions(i,4) = -constant; % tx (pressure)
            boundaryConditions(i,5) = 0; % ty (no shear)
            boundaryConditions(i,6) = 0; % tz (no shear)
        elseif ismember(i, outletNodes)
            % Outlet: pressure and shear 
            boundaryConditions(i,4) = constant; % tx (pressure)
            boundaryConditions(i,5) = 0; % ty (no shear)
            boundaryConditions(i,6) = 0; % tz (no shear)
        elseif ismember(i, [bottomNodes, topNodes])
            % Bottom/Top walls
            total = 0;
            for j = 1:2:9
                total = total + constant * (-1)^((j-1)/2) * (1-(cosh(j*pi*p_z/(2*a))/cosh(j*pi*b/(2*a)))) * (-sin(j*pi*p_y/(2*a))*j*pi/(2*a))/j^3;
            end
            boundaryConditions(i,4) = total*viscosity; % tx (shear)
            boundaryConditions(i,5) = 0; % ty (normal)
            boundaryConditions(i,6) = 0; % tz (shear)
        else
            % Side walls
            boundaryConditions(i,4) = 0; % tx (shear)
            boundaryConditions(i,5) = 0; % ty (shear)
            boundaryConditions(i,6) = 0; % tz (normal)
        end
    end
end

function boundaryConditions = applyMixedBoundaries(boundaryConditions, p, a, b, dp_dx, viscosity, inletNodes, outletNodes, rightNodes, leftNodes, bottomNodes, topNodes)
    constant = 16*a^2/(viscosity*pi^3)*(-dp_dx);
    
    % Apply to all nodes first
    for i = 1:size(p,1)
        boundaryConditions(i,2) = 0; % v=0 everywhere
        boundaryConditions(i,3) = 0; % w=0 everywhere
        
        p_y = p(i,2) - a;
        p_z = p(i,3) - b;
        
        if ismember(i, inletNodes)
            % Inlet: Dirichlet for velocity
            total = 0;
            for j = 1:2:9
                total = total + constant * (-1)^((j-1)/2) * (1-(cosh(j*pi*p_z/(2*a))/cosh(j*pi*b/(2*a)))) * cos(j*pi*p_y/(2*a))/j^3;
            end
            boundaryConditions(i,1) = total;

        elseif ismember(i , outletNodes)
            total = 0;
            totalOutlet = 0;
            totalZ = 0;
            boundaryConditions(i,1) = total;
            for j = 1:2:9
                temp = (-1)^((j-1)/2);
                % ux
                totalOutlet = totalOutlet + constant * temp * (1-(cosh(j*pi*p_z/(2*a))/cosh(j*pi*b/(2*a)))) * cos(j*pi*p_y/(2*a))/j^3;
                
                % ty = du/dy
                total = total + constant * temp * (1-(cosh(j*pi*p_z/(2*a))/cosh(j*pi*b/(2*a)))) * (-sin(j*pi*p_y/(2*a))*j*pi/(2*a))/j^2;
               
                % tz = du/dz
                totalZ = totalZ + constant * temp * (-sinh(j*pi*p_z/(2*a))/cosh(j*pi*b/(2*a)))*j*pi/(2*a) * cos(j*pi*p_y/(2*a))/j^3;
            end
            boundaryConditions(i,1) = totalOutlet; % -p_x value not known
            boundaryConditions(i,2:4) = NaN;
            boundaryConditions(i,5) = total;
            boundaryConditions(i,6) = totalZ;

        elseif ismember(i, [bottomNodes, topNodes])
            % Bottom/Top: no-slip Dirichlet
            boundaryConditions(i,1) = 0;
        elseif ismember(i, [leftNodes, rightNodes])
            % Side walls: no-slip Dirichlet
            boundaryConditions(i,1) = 0;
        end
    end
end