function points = rotatePoints (points, midPoint, angle)
    
    % Init
    n = size(points, 1);
    
    % Shift points by -midpoint
    temp = points - ones(n,1)*midPoint;
    
    % Rotate
    s = sin(angle);
    c = cos(angle);
    temp = [ temp(:,1)*c-temp(:,2)*s, temp(:,2)*c+temp(:,1)*s ];
    
    % Shift back by midpoint
    points = temp + ones(n,1)*midPoint;

end