% A : quantity to be interpolated in
% x : first component, generally the density
% y : second component, generally the specific internal energy
% X : Grid vector for the first component
% Y : Grid vector for the second component

function value = BilinearInterpolation(A,x,y,X,Y)

    % First, interpolate into grid vectors D and U to find the four
    % adjacent points to our query point (x,y)
    DistX = abs(X - x); % component-wise
    DistY = abs(Y - y);
    [~,Ix] = min(DistX);
    
    x1 = X(Ix);
    x2 = X(Ix + 1);
    Idx1 = Ix;
    Idx2 = Ix + 1;
    
    if x - X(Ix) < 0
        x1 = X(Ix - 1);
        x2 = X(Ix);
        Idx1 = Ix - 1;
        Idx2 = Ix;
    end
    
    [~,Iy] = min(DistY);

    
    y1 = Y(Iy);
    y2 = Y(Iy + 1);
    Idy1 = Iy;
    Idy2 = Iy + 1;
    
    if y - Y(Iy) < 0
        y1 = Y(Iy - 1);
        y2 = Y(Iy);
        Idy1 = Iy - 1;
        Idy2 = Iy;
    end
        
    Dx = x2 - x1;
    Dy = y2 - y1;
    
    dx = x - x1;
    dy = y - y1;
    
    Dfx = A(Idx2,Idy1) - A(Idx1,Idy1);
    Dfy = A(Idx1,Idy2) - A(Idx1,Idy1);
    Dfxy = A(Idx1,Idy1) + A(Idx2,Idy2) - A(Idx2,Idy1) - A(Idx1,Idy2);
    
    value = Dfx*dx/Dx + Dfy*dy/Dy + Dfxy*dx/Dx*dy/Dy + A(Idx1,Idy1);
    
end