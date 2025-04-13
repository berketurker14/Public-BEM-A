function plotGeometryAndConnectivity(p, t)

figure;
hold on;
axis equal;
view(3);
xlabel('X'); ylabel('Y'); zlabel('Z');
title('3D Surface Mesh Connectivity');

for elem = 1:size(t,1)
    nodes = t(elem,:);
    X = p(nodes,1);
    Y = p(nodes,2);
    Z = p(nodes,3);
    patch(X,Y,Z,'cyan','FaceAlpha',0.5,'EdgeColor','black');
end

plot3(p(:,1),p(:,2),p(:,3),'ro','MarkerSize',1,'MarkerFaceColor','r');

end
