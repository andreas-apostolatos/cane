classdef testIGAComputationalFluidDynamicsAnalysis < matlab.unittest.TestCase
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Class definition
%
% Test suites for the isogeometric incompressible flow analysis.
%
%% Method definitions
methods (Test)
    testIGA4StokesSteadyState2D(testCase)
    testIGA4TransientTaylorGreenVortices2D(testCase)
end

end
