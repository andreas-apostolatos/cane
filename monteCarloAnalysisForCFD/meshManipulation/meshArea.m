function area = meshArea(nodes, elements)
    % A = meshArea( nodes, elements )
    % Sums up the areas of all elements in the mesh

    % Init
    area = 0;
    
    % Add up all the element areas
    for k=1:length(elements)
        area = area + triangleArea( nodes(elements(k,1:3),1:2) );
    end

end