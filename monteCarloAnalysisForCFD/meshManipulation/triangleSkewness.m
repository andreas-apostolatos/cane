function skewness = triangleSkewness( points )
    % skewness = triangleSkewness(points)
    % Computes the equilateral-area-based skewness of the element defined 
    % by the input 2D points
    % skewness = ( equilateral_area - area ) / equilateral_area
    
    % Init
    if size(points,2) > size(points,1)
        points = points';
    end
    
    % Get triangle area
    area = triangleArea(points);
    
    % Get the circumradius from area ( r = abc/4A )
    radius =                                                            ...
        prod( vecnorm( diff([points;points(1,:)])' ) )                  ...
        /                                                               ...
        ( 4*area );
    
    % Get equilateral area
    equilateralArea = radius*radius * sqrt(3)*3/4;
    
    % Compute skewness
    skewness = (equilateralArea - area) / equilateralArea;

end