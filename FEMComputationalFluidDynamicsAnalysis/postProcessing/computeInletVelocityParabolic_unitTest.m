function valuesInhomDBCModified = ...
    computeInletVelocityParabolic_unitTest ...
    (fldMsh, inhomDBC, valuesInhomDBC, Umid)
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
% Returns modified inhomogeneous boundary conditions, works only for 2D 
% rectangular meshes 
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

inletNodes = floor(inhomDBC/3) + mod(inhomDBC,3 );

% Sort according to y-coordinate
[sortedCoordinates, rowIndices] = sort(fldMsh.nodes(inletNodes,2),'ascend');

% Find the inlet height and relative y-coordinate of each node
H = sortedCoordinates(end) - sortedCoordinates(1);
y = sortedCoordinates - sortedCoordinates(1); 

% Calculate input velocities based on parabolic distribution formula
values = 4*Umid*y.*(H-y)/(H*H);

% Insert values
valuesInhomDBCModified = valuesInhomDBC;
valuesInhomDBCModified(rowIndices) = values;
    
end
