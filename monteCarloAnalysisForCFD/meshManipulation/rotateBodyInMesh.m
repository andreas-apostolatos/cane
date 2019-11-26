function newNodes = rotateBodyInMesh ( mesh, bodyNodeIndices, angle, angleExponent, numberOfNodeLayers )
    % rotatedNodes = rotateBodyInMesh( problemStruct, bodyNodeIndices, angle, exponent, numberOfNodeLayers )
    % Rotates the nodes of bodyNodeIndices by angle, and all connected
    % nodes up to numberOfNodeLayers connection level, by an exponentially
    % decreasing angle. The center of rotation is the mean of the body
    % nodes.
    % The exponent argument should be a 2x1 vector specifying the
    % maximum and minimum exponents respectively (the exponents will 
    % decrease in a quadratic trend, with 0 initial slope)
    %
    % This function performs automatic mesh overlap checks. If the mesh has
    % overlaps, try playing around with different numberOfNodeLayers levels and
    % exponent bounds, and check how the mesh changes.

    % Get initial mesh area (check for overlapping at the end)
    Amsh = meshArea( mesh.nodes, mesh.elements );
    
    % Prealloc
    newNodes = mesh.nodes;

    % Find nodes attached to the body nodes (indirectly if numberOfNodeLayers > 1)
    nodeSets = attachedNodes( bodyNodeIndices, mesh.elements, numberOfNodeLayers );
    
    % Find center of the body
    center = mean( mesh.nodes(bodyNodeIndices, 1:2) );
    
    % Rotate body nodes by the full angle
    newNodes(bodyNodeIndices, 1:2) =                                    ...
        rotatePoints(                                                   ...
            mesh.nodes(bodyNodeIndices, 1:2),                           ...
            center,                                                     ...
            angle                                                       ...
        );
    
    % Quadratically decreasing exponents
    if numberOfNodeLayers ~= 1
        A = [ 0 1 2; 1 1 1; 1 numberOfNodeLayers numberOfNodeLayers*numberOfNodeLayers ];
        coefs = ( A\[ 0; angleExponent(1); angleExponent(2) ] )';
        angleExponent = coefs * [ ones(1,numberOfNodeLayers); 1:numberOfNodeLayers; (1:numberOfNodeLayers).^2 ];
    end
    
    % Rotate attached nodes by exponentially diminishing angles
    for k=1:numberOfNodeLayers
        
        % Set new rotation angle
        angle = angle * angleExponent(k);
        % Rotate current set of nodes
        newNodes(nodeSets{k}, 1:2) =                                    ...
            rotatePoints(                                               ...
                mesh.nodes(nodeSets{k}, 1:2),                           ...
                center,                                                 ...
                angle                                                   ...
            );
        
    end
    
    % Check for overlapping
    Amsh_rotated = meshArea( newNodes, mesh.elements );
    if abs(Amsh - Amsh_rotated) > 1e-10
        error( [                                                        ...
            'Error in rotateBodyInMesh: rotated mesh overlapping ',     ...
            num2str( (Amsh_rotated-Amsh)/Amsh*100 )                     ...
            '%'                                                         ...
            ] )
    end
    
    % Compute skewness
    skewness = meshSkewness( newNodes, mesh.elements );
    if skewness > 0.75
        warning( [                                                      ...
            'Poor mesh quality after rotation! Mesh skewness: ',        ...
            num2str(skewness)                                           ...
            ] );
    end
end