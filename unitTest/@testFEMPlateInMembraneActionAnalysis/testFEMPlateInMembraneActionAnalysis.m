%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
classdef testFEMPlateInMembraneActionAnalysis < matlab.unittest.TestCase
    %% Class definition
    %
    % Test suites for the nonlinear steady-state and transient finite
    % element formulation for the plate in membrane action problem.
    %
    %% Method definitions
    methods (Test)
        testCurvedPlateInMembraneActionSteadyStateLinear(testCase)
        testCantileverBeamTransient(testCase)
    end
end
