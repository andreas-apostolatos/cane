function valuesInhomDBCModified = ...
    computeInletVelocityPowerLaw(fldMsh, inhomDBC, valuesInhomDBC, u0)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the velocity at the nodes on the inlet boundary based on a power
% law assuming zero velocity at the node with the lowest height and maximum
% velocity at the node with the tallest height.
%
%                  Input :
%                 fldMsh : Nodes and elements for the fluid mesh
%               inhomDBC : The global numbering of the nodes where 
%                          inhomogeneous Dirichlet boundary conditions are 
%                          applied
%         valuesInhomDBC : The prescribed values for the inhomogeneous 
%                          Dirichlet boundary conditions
%                     u0 : Velocity at the node with the tallest height
%
%                 Output :
% valuesInhomDBCModified : Array containing the updated prescribed values 
%                          at the nodes where inhomogeneous Dirichlet 
%                          boundary conditions are applied based on the
%                          power law
%
%% Function main body

n = 1/7;

inletNodes = floor(inhomDBC/3) + mod(inhomDBC,3 );

% Sort according to y-coordinate
[sortedCoordinates, rowIndices] = sort(fldMsh.nodes(inletNodes,2),'ascend');

% Find the inlet height and relative y-coordinate of each node
y0 = sortedCoordinates(1);
H = sortedCoordinates(end) - y0;
y = sortedCoordinates - sortedCoordinates(1); 

% Calculate input velocities based on the power law
values = 4*u0*((y - y0)./H).^n;

% Insert values
valuesInhomDBCModified = valuesInhomDBC;
valuesInhomDBCModified(rowIndices) = values;
    
end