classdef testFEMThermalConductionAnalysis < matlab.unittest.TestCase
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Class definition
%
% Test suites for the linear transient thermal conduction analysis.
%
%% Method definitions
methods (Test)
    testTransientSquareCavity(testCase)
    testTransientWallHeating(testCase)
end

end