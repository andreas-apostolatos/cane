function indices = singularIndices(output)
    % indices = singularIndices(output)
    % Returns the indices of ridiculusly large components in 'output'
    % 'output' must be a vector

    % Init
    if size(output,2)>size(output,1)
        output = output';
    end
    origIndices = (1:length(output))';
    
    % Filter NaNs
    % Check nan values
    NaNs        = isnan(output);
    output      = output(~NaNs);
    origIndices = origIndices(~NaNs);
    
    % Get original mean
    output      = abs(output);
    avg         = mean(output);
    % Weed out the big components and compute the new mean
    indices     = output < 1e16;
    newOutput   = output( indices );
    newAvg      = mean(newOutput);
    % Check if the new average is drastically smaller
    if newAvg/avg<1e-5
        indices = find(~indices);
        indices = origIndices(indices);
    else
        indices = [];
    end
    % Add NaNs
    indices = [indices;find(NaNs)];
    

end