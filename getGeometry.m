function [p, t] = getGeometry(Lx, Ly, Lz, Nx, Ny, Nz)

    node_num = 0;
    p = [];
    % Preallocations
    node_num_face1 = zeros(Ny+1, Nz+1); 
    node_num_face2 = zeros(Ny+1, Nz+1); 
    node_num_face3 = zeros(Nx+1, Nz+1);
    node_num_face4 = zeros(Nx+1, Nz+1);
    node_num_face5 = zeros(Nx+1, Ny+1);
    node_num_face6 = zeros(Nx+1, Ny+1);
    
    delta_x = Lx / Nx;
    delta_y = Ly / Ny;
    delta_z = Lz / Nz;
    
    % Face 1 (x=0)
    for i = 0:Ny
        for j = 0:Nz
            node_num = node_num + 1;
            x = 0;
            y = i * delta_y;
            z = j * delta_z;
            p(node_num,:) = [x,y,z];
            node_num_face1(i+1,j+1) = node_num;
        end
    end
    
    % Face 2 (x=Lx)
    for i=0:Ny
        for j=0:Nz
            node_num = node_num +1;
            x = Lx; y = i*delta_y; z = j*delta_z;
            p(node_num,:)=[x,y,z];
            node_num_face2(i+1,j+1)=node_num;
        end
    end
    
    % Face 3 (y=0)
    for i=0:Nx
        for j=0:Nz
            node_num=node_num+1;
            x=i*delta_x; y=0; z=j*delta_z;
            p(node_num,:)=[x,y,z];
            node_num_face3(i+1,j+1)=node_num;
        end
    end
    
    % Face 4 (y=Ly)
    for i=0:Nx
        for j=0:Nz
            node_num=node_num+1;
            x=i*delta_x;y=Ly;z=j*delta_z;
            p(node_num,:)=[x,y,z];
            node_num_face4(i+1,j+1)=node_num;
        end
    end
    
    % Face 5 (z=0)
    for i=0:Nx
        for j=0:Ny
            node_num=node_num+1;
            x=i*delta_x;y=j*delta_y;z=0;
            p(node_num,:)=[x,y,z];
            node_num_face5(i+1,j+1)=node_num;
        end
    end
    
    % Face 6 (z=Lz)
    for i=0:Nx
        for j=0:Ny
            node_num=node_num+1;
            x=i*delta_x;y=j*delta_y;z=Lz;
            p(node_num,:)=[x,y,z];
            node_num_face6(i+1,j+1)=node_num;
        end
    end
    
    % Remove duplicates and orders correctly
    [unique_p, ~, IC] = unique(p, 'rows', 'stable');
    p = unique_p;
    node_num_face1 = IC(node_num_face1);
    node_num_face2 = IC(node_num_face2);
    node_num_face3 = IC(node_num_face3);
    node_num_face4 = IC(node_num_face4);
    node_num_face5 = IC(node_num_face5);
    node_num_face6 = IC(node_num_face6);

    t = [];
    element_num=0;
    
    % Face 1 (Inlet)
    for i=1:Ny
        for j=1:Nz
            n1=node_num_face1(i,j);
            n2=node_num_face1(i,j+1);
            n3=node_num_face1(i+1,j+1);
            n4=node_num_face1(i+1,j);
            element_num=element_num+1;
            t(element_num,:)=[n1,n2,n3,n4];
        end
    end
    
    % Face 2 (Outlet)
    for i=1:Ny
        for j=1:Nz
            n1=node_num_face2(i,j);
            n2=node_num_face2(i,j+1);
            n3=node_num_face2(i+1,j+1);
            n4=node_num_face2(i+1,j);
            element_num=element_num+1;
            t(element_num,:)=[n1,n2,n3,n4];
        end
    end
    
    % Face 3 (Bottom)
    for i=1:Nx
        for j=1:Nz
            n1=node_num_face3(i,j);
            n2=node_num_face3(i,j+1);
            n3=node_num_face3(i+1,j+1);
            n4=node_num_face3(i+1,j);
            element_num=element_num+1;
            t(element_num,:)=[n1,n2,n3,n4];
        end
    end
    
    % Face 4 (Top)
    for i=1:Nx
        for j=1:Nz
            n1=node_num_face4(i,j);
            n2=node_num_face4(i,j+1);
            n3=node_num_face4(i+1,j+1);
            n4=node_num_face4(i+1,j);
            element_num=element_num+1;
            t(element_num,:)=[n1,n2,n3,n4];
        end
    end
    
    % Face 5 (Left)
    for i=1:Nx
        for j=1:Ny
            n1=node_num_face5(i,j);
            n2=node_num_face5(i,j+1);
            n3=node_num_face5(i+1,j+1);
            n4=node_num_face5(i+1,j);
            element_num=element_num+1;
            t(element_num,:)=[n1,n2,n3,n4];
        end
    end
    
    % Face 6 (Right)
    for i=1:Nx
        for j=1:Ny
            n1=node_num_face6(i,j);
            n2=node_num_face6(i,j+1);
            n3=node_num_face6(i+1,j+1);
            n4=node_num_face6(i+1,j);
            element_num=element_num+1;
            t(element_num,:)=[n1,n2,n3,n4];
        end
    end

end
