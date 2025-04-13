function desiredDir = getDesiredDirection(x1, x2, x3, x4, Lx, Ly, Lz) % Not done. Learn further methods. This one's from gpt.

avgX = mean([x1(1), x2(1), x3(1), x4(1)]);
avgY = mean([x1(2), x2(2), x3(2), x4(2)]);
avgZ = mean([x1(3), x2(3), x3(3), x4(3)]);

tol = 1e-6;

if abs(avgX) < tol
    % Inlet face
    desiredDir = [-1, 0, 0];
elseif abs(avgX - Lx) < tol
    % Outlet face
    desiredDir = [1, 0, 0];
elseif abs(avgY - 0) < tol
    % Bottom face
    desiredDir = [0, -1, 0];
elseif abs(avgY - Ly) < tol
    % Top face
    desiredDir = [0, 1, 0];
elseif abs(avgZ - 0) < tol
    % Left face
    desiredDir = [0, 0, -1];
elseif abs(avgZ - Lz) < tol
    % Right face
    desiredDir = [0, 0, 1];
else

    desiredDir = [1,0,0]; % Default or handle differently
end
end
