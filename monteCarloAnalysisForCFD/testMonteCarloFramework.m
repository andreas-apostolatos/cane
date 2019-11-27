function testMonteCarloFramework(testCase)
    
    %% INIT
    caseName    = 'first_test';
    tolerance   = 1e-15;

    %% SYSTEM TEST
    % Define input and output functions
    inputFunction   = @input_inletVelocity;
    outputFunction  = @(paramStruct,up,FComplete) sum(abs(FComplete));
    
    % Define input values
    input           = linspace(0,1,5);
    
    % Run calculation
    [output, problemStruct] = monteCarlo(                               ...
        caseName,                                                       ...
        input,                                                          ...
        inputFunction,                                                  ...
        outputFunction                                                  ...
        );
    
    controlValues = [                                                   ...
        0;                                                              ...
        4.26036030572575;                                               ...
        16.6134945284118;                                               ...
        37.0723416125103;                                               ...
        65.6359940356688];
    
    testCase.verifyEqual(                                               ...
        cell2mat(output), controlValues,                                ...
        'AbsTol',tolerance);
    
    
end