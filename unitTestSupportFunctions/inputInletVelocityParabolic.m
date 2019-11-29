%% Function documentation
%
% Returns modified boundary conditions, works only for 2D rectangular meshes 
% Example: we prescbribe input velocity on left boundary, function changes that
% velocity to have parabolic distribution with new Umax
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
function valuesInhomDBCModified = ...
    inputInletVelocityParabolic(fldMsh, inhomDBC, valuesInhomDBC, Umid)
   
    % Find inlet
    rowIndices = findInletNodeIndices(fldMsh, inhomDBC);
    rowIndices = inhomDBC( rowIndices );
    rowIndices = floor(rowIndices/3) + mod(rowIndices,3);
    
    % Sort 
    [sortedCoordinates, rowIndices] = sort(                             ...
        fldMsh.nodes(rowIndices,2),'ascend');
    
    % Get parabolic distribution
    H = sortedCoordinates(end)-sortedCoordinates(1);
    values = parabolicInlet(                                            ...
        sortedCoordinates-sortedCoordinates(1),H,Umid);
    
    % Insert
    valuesInhomDBCModified = valuesInhomDBC;
    valuesInhomDBCModified(rowIndices) = values;
   

    %% nested function that finds indicied of nodes on the inlet
    function indices = findInletNodeIndices(fldMsh, inhomDBC)
        rows = floor(inhomDBC/3) + mod(inhomDBC,3 );
        indices = find(                                                 ... 
            abs(fldMsh.nodes(rows,1) - min(fldMsh.nodes(:,1))) < 1e-10);
    end
    
    %% nested function that defines the parabolic distribution
    function value = parabolicInlet(y,H,Um)
        value = 4*Um*y.*(H-y)/(H*H);
    end
    
end