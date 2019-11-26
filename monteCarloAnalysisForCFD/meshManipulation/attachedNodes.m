function nodes = attachedNodes (baseNodes, elements, numberOfNodeLayers)

    % Prealloc
    nodes = cell(numberOfNodeLayers, 1);
    
    % Find indices connected to the base nodes, one by one per node
    for k1=1:numberOfNodeLayers
        for k2=1:length(baseNodes)
            
            % Find continous indices
            indices = find(elements == baseNodes(k2));
            % Convert continous indices to element numbers
            indices = mod( indices, size(elements,1) );
            % Check if element was used before
            if elements(indices, 1) ~= 0
                % Add new nodes to the list
                for k3=1:length(indices)
                    newNodes = elements(indices(k3),:);
                    newNodes = newNodes( newNodes ~= baseNodes(k2) );
                    nodes{k1} = [ nodes{k1}, newNodes ];
                end
                % Set used elements to 0
                elements(indices, :) = 0;
            end
            
        end
        
        % Delete base nodes from the list
        for k2=1:length(baseNodes)
            nodes{k1} = nodes{k1}( nodes{k1} ~= baseNodes(k2) );
        end
        
        % Delete repeated nodes
        nodes{k1} = unique(nodes{k1});
        
        % Update baseNodes
        baseNodes = nodes{k1};
        
    end

end