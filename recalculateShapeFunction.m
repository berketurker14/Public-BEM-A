function N=recalculateShapeFunction(gp)
    % Linear 2D For now
    N = [0.5 * (1 - gp) , 0.5 * (1 + gp)];
end