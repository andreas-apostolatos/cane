function Fy = lift (paramStruct, up, FComplete)

    F = bodyForces(paramStruct, FComplete);
    Fy = sum(F(:,2));
    
end