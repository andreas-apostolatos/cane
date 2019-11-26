function [] = input_viscosity( paramStruct, input )
    % [] = input_viscosity( paramStruct, input )
    % Sets the dynamic viscosity of the paramStruct to input
    
    paramStruct.setField(                                               ...
        {'parameters.nue'},                                             ...
        {input}                                                         ...
        );


end