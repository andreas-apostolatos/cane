%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
classdef testIGAKirchhoffLoveShellAnalysis < matlab.unittest.TestCase
    %% Class definition
    %
    % Test suites for the linear, nonlinear and multipatch isogeometric 
    % Kirchhoff-Love shell analysis.
    
    %% Method definitions
    methods (Test)
        testLinearKirchoffLoveShellAnalysis(testCase)
        testNonlinearKirchoffLoveShellAnalysis(testCase)
        testPlateBentIntoACircleKirchoffLoveShellAnalysis(testCase)
        testLinearKirchoffLoveShellMultipatchAnalysis(testCase)
        testNonlinearKirchoffLoveShellMultipatchAnalysis(testCase)
    end
end
