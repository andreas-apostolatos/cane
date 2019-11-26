function [] = input_density( paramStruct, input )
    % [] = input_density( paramStruct, input )
    % Sets the density of the paramStruct to input
    
    paramStruct.setField(                                               ...
        {'parameters.rho'},                                             ...
        {input}                                                         ...
        );

end