function problemStruct = init(caseName, pathToCase)
    % problemStruct = init(varargin)
    % Changes the working directory to root/monteCarlo and returns a
    % problemStructure (structure that contains all parameters/settings 
    % needed to run an analysis).
    % Input argument:    - name of the case .dat file that gets loaded
    %                    - path to the case where the file is   
    % Case files are loaded from:
    % root/inputGiD/FEMComputationalFluidDynamicsAnalysis

    %% Includes
    % Get the full path to this function
    rootFolder  = mfilename('fullpath');
    
    % Get the path to the root folder (1 level up)
    rootFolder(rootFolder=='\') = '/';
    folderIndices   = strfind(rootFolder,'/');
    rootFolder      = rootFolder(1:folderIndices(end-1)-1);
    
    % Add the Monte Carlo folders
    addpath([rootFolder,'/monteCarloAnalysisForCFD'])
    addpath([rootFolder,'/monteCarloAnalysisForCFD/inputFunctions'])
    addpath([rootFolder,'/monteCarloAnalysisForCFD/postprocessing'])
    addpath([rootFolder,'/monteCarloAnalysisForCFD/meshManipulation'])
        
    % Add functions related to equation system solvers
    addpath([rootFolder,'/equationSystemSolvers/']);

    % Add general math functions
    addpath([rootFolder,'/generalMath/']);

    % Add the classical finite element basis functions
    addpath([rootFolder,'/basisFunctions/']);

    % Add all functions related to the Finite Element Methods for Computational
    % Fluid Dynamics problems
    addpath([rootFolder,'/FEMComputationalFluidDynamicsAnalysis/solutionMatricesAndVectors/'],...
            [rootFolder,'/FEMComputationalFluidDynamicsAnalysis/initialConditions'],...
            [rootFolder,'/FEMComputationalFluidDynamicsAnalysis/solvers/'],...
            [rootFolder,'/FEMComputationalFluidDynamicsAnalysis/loads/'],...
            [rootFolder,'/FEMComputationalFluidDynamicsAnalysis/output/'],...
            [rootFolder,'/ALEMotion/']);

    % Add functions for the mesh motion 
    addpath([rootFolder,'/FEMPlateInMembraneActionAnalysis/solutionMatricesAndVectors/']);

    % Add all functions related to parsing
    addpath([rootFolder,'/parsers/']);

    % Add all functions related to the efficient computation functions
    addpath([rootFolder,'/efficientComputation/']);

    % Add all functions related to ALE
    addpath([rootFolder,'/ALEMotion/']);

    %% GUI
    % On the body forces
    problemStruct.computeBodyForces = @computeConstantVerticalBodyForceVct;

    % On the initial conditions
    % computeInitialConditions = @computeInitialConditionsFromVTKFileFEM4NSE2D;
    problemStruct.computeInitialConditions = @computeNullInitialConditionsFEM4NSE2D;

    % Case name
    problemStruct.pathToCase = pathToCase;
    problemStruct.caseName = caseName;
  

    %% Parse the data from the GiD input file
    [
        problemStruct.fldMsh,                                           ...
        problemStruct.homDBC,                                           ...
        problemStruct.inhomDBC,                                         ...
        problemStruct.valuesInhomDBC,                                   ...
        problemStruct.nodesALE,                                         ...
        problemStruct.NBC,                                              ...
        problemStruct.analysis,                                         ...
        problemStruct.parameters,                                       ...
        problemStruct.propNLinearAnalysis,                              ...
        problemStruct.propFldDynamics,                                  ...
        problemStruct.gaussInt                                          ...
    ]                                                                   ...
    = parse_FluidModelFromGid(                                          ...
        problemStruct.pathToCase,                                       ...
        problemStruct.caseName,                                         ...
        'outputEnabled'                                                 ...
    );

    %% Choose the equation system solver
    if strcmp(problemStruct.analysis.type,'NAVIER_STOKES_2D')
        problemStruct.solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;
    elseif strcmp(problemStruct.analysis.type,'NAVIER_STOKES_3D')
        problemStruct.solve_LinearSystem = @solve_LinearSystemGMResWithIncompleteLUPreconditioning;
    else
        error('Neither NAVIER_STOKES_2D or NAVIER_STOKES_3D has been chosen');
    end

    %% Define the name of the vtk file from where to resume the simulation
    problemStruct.VTKResultFile = 'undefined';
end