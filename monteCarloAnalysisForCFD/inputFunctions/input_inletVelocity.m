function [] = input_inletVelocity( paramStruct, input )
    % [] = input_inletVelocity( paramStruct, input )
    % Sets the inlet velocity of the paramStruct to input
    
    paramStruct.setField(                                               ...
        {'valuesInhomDBC'},                                             ...
        {input*ones(size(paramStruct.valuesInhomDBC))});

end