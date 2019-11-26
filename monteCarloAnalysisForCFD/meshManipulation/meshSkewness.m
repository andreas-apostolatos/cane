function [maxSkewness, skewness] = meshSkewness( nodes, elements )
    % [maxSkewness, skewness] = meshSkewness( nodes, elements )
    % Computes the skewness of every element in a triangular mesh

    skewness = zeros( length(elements),1 );
    
    for k=1:length(elements)
        skewness(k) = triangleSkewness( nodes(elements(k,:),1:2) );
    end
    
    maxSkewness = max(skewness);

end