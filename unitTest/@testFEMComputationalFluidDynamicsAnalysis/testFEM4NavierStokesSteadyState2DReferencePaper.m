function testFEM4NavierStokesSteadyState2DReferencePaper(testCase)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Tests the stabilized finite element analysis for the 2D steady-state
% Navier-Stokes equations, over the flow around a cylinder problem with 
% respect to reference paper by:
% http://www.mathematik.tu-dortmund.de/lsiii/cms/papers/SchaeferTurek1996.pdf
%
% Function layout :
%
% 0. Read input
%
% 1. Parse the data from the GiD input file
%
% 2. GUI
%
% 3. Choose the equation system solver
%
% 4. Define the name of the vtk file from where to resume the simulation
%
% 5. Solve the CFD problem
%
% 6. Calculate drag and lift force from the nodal forces
%
% 7. Calculate drag and lift coefficiend based on drag and lift force
%
% 8. Define the expected solutions
%
% 9. Verify the results
%
%% Function main body

%% 0. Read input

% Define relative tolerance
relTol = 0.01;

% Define the path to the case
pathToCase = '../../inputGiD/FEMComputationalFluidDynamicsAnalysis/';
caseName = 'unitTest_flowAroundCylinder2DReferencePaper';

%% 1. Parse the data from the GiD input file
problemStruct = init(caseName, pathToCase);
paramStruct = ParameterStructure(problemStruct);

% [fldMsh,homDBC,inhomDBC,valuesInhomDBC,nodesALE,~,analysis,parameters,...
%     propNLinearAnalysis,propFldDynamics,gaussInt] = ...
%     parse_FluidModelFromGid...
%     (pathToCase,caseName,'');

%% 2. GUI

% On the body forces
%computeBodyForces = @computeConstantVerticalBodyForceVct;

% On the initial conditions
% computeInitialConditions = @computeInitialConditionsFromVTKFileFEM4NSE2D;
%computeInitialConditions = @computeNullInitialConditionsFEM4NSE2D;


% % On the transient analysis properties
% if strcmp(propFldDynamics.method,'bossak')
%     propFldDynamics.computeProblemMtrcsTransient = ...
%         @computeProblemMtrcsBossakFEM4NSE;
%     propFldDynamics.computeUpdatedVct = ...
%         @computeBossakTIUpdatedVctAccelerationFieldFEM4NSE2D;
% end

%% 3. Choose the equation system solver
% if strcmp(analysis.type,'NAVIER_STOKES_2D')
%     solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;
% elseif strcmp(analysis.type,'NAVIER_STOKES_3D')
%     solve_LinearSystem = @solve_LinearSystemGMResWithIncompleteLUPreconditioning;
% else
%     error('Neither NAVIER_STOKES_2D or NAVIER_STOKES_3D has been chosen');
% end

%% 4. Define the name of the vtk file from where to resume the simulation
% VTKResultFile = 'undefined';


%% 4. Change the input velocity to match the paper - parabolic input

% max input velocity defined in the reference paper
Umax = 0.3;
% change the input velocity to have the parabolic distribution
input_inletVelocityParabolic(paramStruct, Umax);


%% 5. Solve the CFD problem
% [up,FComplete,hasConverged,~] = solve_FEMVMSStabSteadyStateNSE2D...
%     (fldMsh,homDBC,inhomDBC,valuesInhomDBC,nodesALE,parameters,...
%     computeBodyForces,analysis,computeInitialConditions,...
%     VTKResultFile,solve_LinearSystem,propFldDynamics,propNLinearAnalysis,...
%     gaussInt,caseName,'');

fieldSelect = @fieldSelector;
[up,FComplete,hasConverged, ~] = ...
            solve_FEMVMSStabSteadyStateNSE2D(                                          ...
            feval(fieldSelect,paramStruct,problemStruct,'fldMsh'),                   ...
            feval(fieldSelect,paramStruct,problemStruct,'homDBC'),                   ...
            feval(fieldSelect,paramStruct,problemStruct,'inhomDBC'),                 ...
            feval(fieldSelect,paramStruct,problemStruct,'valuesInhomDBC'),           ...
            feval(fieldSelect,paramStruct,problemStruct,'nodesALE'),                 ...
            feval(fieldSelect,paramStruct,problemStruct,'parameters'),               ...
            feval(fieldSelect,paramStruct,problemStruct,'computeBodyForces'),        ...
            feval(fieldSelect,paramStruct,problemStruct,'analysis'),                 ...
            feval(fieldSelect,paramStruct,problemStruct,'computeInitialConditions'), ...
            feval(fieldSelect,paramStruct,problemStruct,'VTKResultFile'),            ...
            feval(fieldSelect,paramStruct,problemStruct,'solve_LinearSystem'),       ...
            feval(fieldSelect,paramStruct,problemStruct,'propFldDynamics'),          ...
            feval(fieldSelect,paramStruct,problemStruct,'propNLinearAnalysis'),      ...
            feval(fieldSelect,paramStruct,problemStruct,'gaussInt'),                 ...
            feval(fieldSelect,paramStruct,problemStruct,'caseName'),                 ...
            'outputDisabled'                                                            ...
        );

%% 6. Calculate drag and lift force from the nodal forces

Fx = drag (paramStruct, up, FComplete);
Fy = lift (paramStruct, up, FComplete);


%% 7. Calculate drag and lift coefficient based on drag and lift force

% define parameters used in reference paper and simualiton
ro = 1;
Ubar = 0.2;
D = 0.1;
% calculate drag and lift coefficiet
dragCoefficient = (2 * Fx )/(ro * Ubar * Ubar * D);
liftCoefficient = (2 * Fy )/(ro * Ubar * Ubar * D);

%% 8. Define the expected solutions

% Define the expected solution in terms of calculated drag coefficient
% 5.5700 - 5.5900 (upper and lower bound in the paper, 2D case)
expSolDragCoefficient = 5.58;
 
% Define the expected solution in terms of calculated lift coefficient
% 0.0104 - 0.0110 (upper and lower bound in the paper, 2D case)
expSolLiftCoefficient = 0.0107;
   
% Define the expected solution in terms of the convergence flag
expSolHasConverged = true;

%% 9. Verify the results
testCase.verifyEqual(dragCoefficient,expSolDragCoefficient,'RelTol',relTol);
testCase.verifyEqual(liftCoefficient,expSolLiftCoefficient,'RelTol',relTol);
testCase.verifyEqual(hasConverged,expSolHasConverged,'RelTol',relTol);

%% FUNCTION DEFINITIONS
    function paramField = fieldSelector(parameterStructure,problemStructure,fieldName)
        % If parameterStructure contains 'fieldName', return it. Otherwise,
        % return 'fieldName' from the problemStructure
        % This function ensures that you don't have to modify the arguments
        % passed to the FE solver when you add/remove fields to/from the 
        % ParameterStructure class
        if isprop(parameterStructure,fieldName)
            paramField = parameterStructure.(fieldName);
        else
            paramField = problemStructure.(fieldName);
        end
    end
end
