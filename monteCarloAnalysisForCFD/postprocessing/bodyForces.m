function F = bodyForces (paramStruct, FComplete)
    
    % Get node indices on the body
    indices = getBodyNodes(paramStruct.fldMsh, paramStruct.homDBC);
    
    % Reshape the input forces to [Fx, Fy, Fz]
    FComplete = -reshape(FComplete, [3, length(FComplete)/3])';
    
    % Select forces on body nodes
    F = FComplete(indices, 1:2) * paramStruct.parameters.rho;

end