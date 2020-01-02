%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _______________________________________________________               %
%   _______________________________________________________               %
%                                                                         %
%   Technische Universität München                                        %
%   Lehrstuhl für Statik, Prof. Dr.-Ing. Kai-Uwe Bletzinger               %
%   _______________________________________________________               %
%   _______________________________________________________               %
%                                                                         %
%                                                                         %
%   Authors                                                               %
%   _______________________________________________________________       %
%                                                                         %
%   Dipl.-Math. Andreas Apostolatos    (andreas.apostolatos@tum.de)       %
%   Altug Emiroglu, M.Sc.              (altug.emiroglu@tum.de)            %
%   Dr.-Ing. Roland Wüchner            (wuechner@tum.de)                  %
%   Prof. Dr.-Ing. Kai-Uwe Bletzinger  (kub@tum.de)                       %
%   _______________________________________________________________       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Script documentation
%
% Task : Geometrically nonlinear transient plane stress analysis
%
% Date : 19.02.2014
%
%% Preamble

% Clear memory
clear;

% Clear the command window
clc;

% Close all generated windows
close all;

%% Includes

% Add general math functions
addpath('../../generalMath/');

% Add all functions related to parsing
addpath('../../parsers/');

% Add all functions related to the low order basis functions
addpath('../../basisFunctions/');

% Equation system solvers
addpath('../../equationSystemSolvers/');

% Equation system solvers
addpath('../../efficientComputation/');

% Include functions related to transient analysis
addpath('../../transientAnalysis/');

% Add all functions related to plate in membrane action analysis
addpath('../../FEMPlateInMembraneActionAnalysis/solvers/',...
        '../../FEMPlateInMembraneActionAnalysis/solutionMatricesAndVectors/',...
        '../../FEMPlateInMembraneActionAnalysis/loads/',...
        '../../FEMPlateInMembraneActionAnalysis/graphics/',...
        '../../FEMPlateInMembraneActionAnalysis/output/',...
        '../../FEMPlateInMembraneActionAnalysis/initialConditions/',...
        '../../FEMPlateInMembraneActionAnalysis/postprocessing/');

%% Parse data from GiD input file

% Define the path to the case
pathToCase = '../../inputGiD/FEMPlateInMembraneActionAnalysis/';

% caseName = 'curvedPlateTipShearPlaneStress';
caseName = 'cantileverBeamPlaneStressTransientNLinear';

% Not a unit test case
isUnitTest = false;

% Parse the data from the GiD input file
[strMsh,homDBC,inhomDBC,valuesInhomDBC,NBC,analysis,parameters,...
    propNLinearAnalysis,propStrDynamics,gaussInt] = ...
    parse_StructuralModelFromGid(pathToCase,caseName,'outputEnabled');

%% GUI

% Equation system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;

% Choose the matric computation corresponding to the chosen time
% integration scheme
if strcmp(propStrDynamics.method,'EXPLICIT_EULER')
    propStrDynamics.computeProblemMtrcsTransient = ...
        @computeProblemMtrcsExplicitEuler;
    propStrDynamics.computeUpdatedVct = ...
        @computeBETITransientUpdatedVctAccelerationField;
elseif strcmp(propStrDynamics.method,'BOSSAK')
    propStrDynamics.computeProblemMtrcsTransient = ...
        @computeProblemMtrcsBossak;
    propStrDynamics.computeUpdatedVct = ...
        @computeBossakTransientUpdatedVctAccelerationField;
else
    error('Invalid time integration method selected in propStrDynamics.method as %s',propStrDynamics.method);
end

% Choose time integration scheme parameters
if strcmp(propStrDynamics.method,'BOSSAK')
    propStrDynamics.alphaB = -.1; % -.1
    propStrDynamics.betaB = .5; % .5
    propStrDynamics.gammaB = .6; % .6
end

% Initialize graphics index
graph.index = 1;

%% Visualization of the configuration
% graph.index = plot_referenceConfigurationFEMPlateInMembraneAction(strMsh,analysis,F,homDBC,graph,'outputEnabled');

%% Solve the plate in membrane action problem
minElSize = solve_FEMPlateInMembraneActionNLinearTransient...
    (analysis,strMsh,homDBC,inhomDBC,valuesInhomDBC,NBC,...
    @computeLoadVctFEMPlateInMembraneAction,...
    parameters,propNLinearAnalysis,propStrDynamics,solve_LinearSystem,...
    gaussInt,caseName,isUnitTest,'outputEnabled');

%% Postprocessing
% graph.visualization.geometry = 'reference_and_current';
% resultant = 'stress';
% component = 'y';
% graph.index = plot_currentConfigurationAndResultants(strMsh,homDBC,dHat,parameters,analysis,resultant,component,graph);

%% END OF THE SCRIPT