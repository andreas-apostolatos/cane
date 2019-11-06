%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
classdef testIGAMembraneAnalysis < matlab.unittest.TestCase
    %% Class definition
    %
    % Test suites for the linear, nonlinear and multipatch isogeometric 
    % membrane analysis.
    
    %% Method definitions
    methods (Test)
        testPagewiseComputationTangStiffMtxIGAMembrane(testCase)
        testFoFiMembraneAnalysis(testCase)
        testNitscheWeakDBCPlateMembraneAnalysis(testCase)
        testNonlinearMembraneAnalysis(testCase)
        testNonlinearMembraneMultipatchAnalysis(testCase)
        testDDMNitschePlateMembraneAnalysis(testCase)
        testTransientNonlinearMembraneAnalysis(testCase)
        testTransientNitscheWeakDBCPlateMembraneAnalysis(testCase)
        testDDMFourPointSailAnalysis(testCase)
    end
end
