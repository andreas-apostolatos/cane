function valuesInhomDBCModified = ...
    computeInletVelocityParabolic_unitTest(fldMsh, inhomDBC, valuesInhomDBC, Umid)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Marko Leskovar
%                  Andreas Apostolatos
%
%% Function documentation
%
% Returns modified inhomogeneous boundary conditions, works only for 2D rectangular meshes 
% Example: we prescbribe input velocity on left boundary, function changes 
% that velocity to have parabolic distribution with new Umid
%
%                  Input :
%                 fldMsh : Nodes and elements for the fluid mesh
%               inhomDBC : The global numbering of the nodes where 
%                          inhomogeneous Dirichlet boundary conditions are 
%                          applied
%         valuesInhomDBC : The prescribed values for the inhomogeneous 
%                          Dirichlet boundary conditions
%                   Umid : Mean value of the velocity of parabolic distribution
%
%                 Output :
% valuesInhomDBCModified : Nodes and elements from the fluid mesh that lie on
%                          a body 
%       
%
%% Function main body

    % find indices of nodes on the inlet
    rows = floor(inhomDBC/3) + mod(inhomDBC,3 );
    rowIndices = abs(fldMsh.nodes(rows,1) - min(fldMsh.nodes(:,1)))     ...
                    < 1e-10;
    
    rowIndices = inhomDBC(rowIndices);
    rowIndices = floor(rowIndices/3) + mod(rowIndices,3);
    
    % Sort 
    [sortedCoordinates, rowIndices] = sort(                             ...
        fldMsh.nodes(rowIndices,2),'ascend');
    
    % calculate input velocities based on parabolic distribution
    H = sortedCoordinates(end)-sortedCoordinates(1);
    y = sortedCoordinates-sortedCoordinates(1); 
    values = 4*Umid*y.*(H-y)/(H*H);
    
    % Insert
    valuesInhomDBCModified = valuesInhomDBC;
    valuesInhomDBCModified(rowIndices) = values;
    
end