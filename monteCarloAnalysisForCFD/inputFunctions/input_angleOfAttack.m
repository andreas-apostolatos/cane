function [] = input_angleOfAttack( paramStruct, input, varargin )
    % [] = input_angleOfAttack( paramStruct, distributionFunction )
    % Changes the <fldMsh.nodes> by rotating the body nodes, and the nodes
    % connected to them up until a numberOfNodeLayers levels. The body
    % nodes are rotated by input.
    % Optional arguments: rotation parameters
    % multiplierRange       : range in which the multipliers can change
    % numberOfNodeLayers    : number of node layers to rotate

    % Rotation parameters
    multiplierRange = [ .99 .55 ];
    numberOfNodeLayers = 12;
    
    % Optional inputs
    if ~isempty(varargin)
        for k=1:numel(varargin)
            if length(varargin{k})>1
                multiplierRange     = varargin{k};
            else
                numberOfNodeLayers  = varargin{2};
            end
        end
    end
       
    % Get body node indices
    bodyNodeIndices = getBodyNodes( paramStruct.fldMsh, paramStruct.homDBC );
    rotatedNodes = rotateBodyInMesh(                                    ...
            paramStruct.fldMsh,                                         ...
            bodyNodeIndices,                                            ...
            input,                                                      ...
            multiplierRange,                                            ...
            numberOfNodeLayers                                          ...
            );
    % Rotate nodes and update structure
    paramStruct.setField(                                               ...
        {'fldMsh.nodes', 'fldMsh.initialNodes'},                        ...
        {rotatedNodes, rotatedNodes}                                    ...
    );

end