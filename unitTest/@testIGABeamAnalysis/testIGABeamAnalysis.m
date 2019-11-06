%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
classdef testIGABeamAnalysis < matlab.unittest.TestCase
    %% Class definition
    %
    % Test suites for the linear Bernoulli and Timoshenko beam analysis
    
    properties
    end
    
    %% Method definitions
    methods (Test)
        testCantileverBeamSubjectToDistributedLoad(testCase)
        testCircularBeamSubjectToInternalPressureLoad(testCase)
    end
end
