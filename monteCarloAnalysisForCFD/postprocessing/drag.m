function Fx = drag (paramStruct, up, FComplete)

    F = bodyForces(paramStruct, FComplete);
    Fx = sum(F(:,1));

end