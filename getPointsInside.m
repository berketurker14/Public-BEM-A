function x_col_inside = getPointsInside(Lx, Ly, Lz, xNum, yNum, zNum, fillPercentage)

    % Fill Percentage -> How close it will be generating near walls 1 -> wall to wall
    
    numInnerPoints = xNum * yNum * zNum;
    x_col_inside = zeros(numInnerPoints,3);
    
    index=1;
    for i=1:xNum
        X = fillPercentage*Lx/xNum*i;
        for j=1:yNum
            Y = fillPercentage*Ly/yNum*j;
            for k=1:zNum
                Z = fillPercentage*Lz/zNum*k;
                x_col_inside(index,:) = [X, Y, Z];
                index=index+1;
            end
        end
    end
end
