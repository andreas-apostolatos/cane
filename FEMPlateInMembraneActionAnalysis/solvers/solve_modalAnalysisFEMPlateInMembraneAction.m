function [eigenmodeShapes, naturalFrequencies, ...
    dHat, FComplete, minElSize] = ...
    solve_modalAnalysisFEMPlateInMembraneAction ...
    (propAnalysis, strMsh, dHat, homDOFs, inhomDOFs, valuesInhomDOFs, ...
    propNBC, computeBodyForces, propParameters, solve_LinearSystem, ...
    propModalAnalysis, propNLinearAnalysis, propGaussInt, propOutput, ...
    caseName, pathToOutput, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the eigenfrequencies and eigenmode shapes field, along with the 
% complete force vector and the minimum element area size corresponding to 
% the modal analysis of a plate in membrane action problem solved with the 
% classical Finite Element Method for the general nonlinear case.
%
%               Input :
%        propAnalysis : Structure containing information about the analysis
%                            .type : The analysis type
%              strMsh : Nodes and elements in the mesh
%                dHat : Initial conditions
%             homDOFs : The global numbering of the nodes where homogeneous
%                       Dirichlet boundary conditions are applied 
%           inhomDOFs : The global numbering of the nodes where 
%                       inhomogeneous Dirichlet boundary conditions are 
%                       applied
%     valuesInhomDOFs : Prescribed values on the nodes where inhomogeneous
%                       Dirichlet boundary conditions are applied
%             propNBC : Structure containing information on the Neumann
%                       boundary conditions
%                          .nodes : The nodes where Neumann boundary 
%                                   conditions are applied
%                       .loadType : The type of the load for each Neumann 
%                                   node
%                      .fctHandle : The function handle for each Neumann 
%                                   node for the computation of the load 
%                                   vector (these functions are unde the 
%                                   folder load)
%   computeBodyForces : Function handle to body force vector computation
%      propParameters : Structure containing the material properties,
%                              .E : Young's modulus
%                            .nue : Poisson's ratio
%                              .t : thickness
%  solve_LinearSystem : Function handle to the solution of the 
%   propModalAnalysis : Structure containing information on the modal
%                       analysis,
%                           .numEig : Number of eigenfrequencies and
%                                     eigenmode shapes to compute
% propNLinearAnalysis : Structure containing information on the
%                       geometrically nonlinear analysis
%                          .method : The employed nonlinear scheme
%                       .tolerance : The residual tolerance
%                         .maxIter : The maximum number of the nonlinear 
%                                    iterations
%       propGaussInt : On the numerical integration (quadrature)
%                       .type : 'default', 'manual'
%                       .noGP : Number of Gauss Points
%         propOutput : Structure containing information on writting the
%                      results for postprocessing,
%                                .isOutput : Flag on whether the results 
%                                            to be written out
%                       .writeOutputToFile : Function handle to the
%                                            writting out of the results
%                           .VTKResultFile : Specifies the name of the
%                                            VTK result file from which
%                                            the simulation to be restarted
%                                            If it is specified as 
%                                            'undefined' the simulation 
%                                            starts from time TStart
%           caseName : The name of the case in the inputGiD case folder
%       pathToOutput : Define the path to where to write out the results
%             outMsg : On outputting information
%
%             Output :
%    eigenmodeShapes : Array containing all the modal shapes
% naturalFrequencies : Vector containing all eigenfrequencies in ascending
%                      order
%               dHat : The nodal displacement field corresponding to the
%                      steady-state solution
%          FComplete : The complete force vector corresponding to the
%                      steady-state solution
%          minElSize : The minimum element area size in the mesh
%
% Function layout :
%
% 0. Read input
%
% 1. Find the prescribed and the free DOFs of the system
%
% 2. Compute the mass matrix of the plate in membrane action
%
% 3. Compute load vector for setting up the initial steady-state computation
%
% 4. Solve the initial steady-state problem
%
% 5. Solve the eigenfrequency problem
%
% 6. Solve the generalized eigenvalue problem to get the eigenmode shapes
%
% 7. Apply the boundary conditions onto the system
%
% 8. Write out the results into a file
%
% 9. Appendix
%
%% Function main body
if ~isfield(propModalAnalysis, 'numEig')
    error('propModalAnalysis must contain field numEig');
else
    if ~isnumeric(propModalAnalysis.numEig)
        error('Number of eigenfrequencies must contain of numeric type');
    else
        if mod(propModalAnalysis.numEig, 1) ~= 0
            error('Number of eigenfrequencies must be an integer');
        else
            if propModalAnalysis.numEig <= 0
                error('Number of eigenfrequencies must be a strictly positive integer');
            end
        end
    end
end
if strcmp(outMsg,'outputEnabled')
    fprintf('__________________________________________________________________\n');
    fprintf('##################################################################\n');
    fprintf('Computation of the eigenfrequencies and mode shapes for a plate in\n');
    fprintf('membrane action has been initiated\n\n');
    fprintf('Number of selected eigenfrequencies : %d \n', propModalAnalysis.numEig);
    fprintf('__________________________________________________________________\n\n');
    tic;
end

%% 0. Read input

% Number of nodes in the mesh
numNodes = length(strMsh.nodes(:, 1));

% Number of DOFs in the mesh
numDOFs = 2*numNodes;

% GLobal DOF numbering
DOFNumbering = 1:numDOFs;

% Assign a dummy load factor
loadFactor = 1;

% Number of eigenfrequencies to be computed
numEig = propModalAnalysis.numEig;

% Assign dummy variables
uSaved = 'undefined';
uDot = 'undefined';
uDDot = 'undefined';
uDotSaved = 'undefined';
uDDotSaved = 'undefined';
propStrDynamics = 'undefined';
uMeshALE = 'undefined';
dDot = 'undefined';
dDDot = 'undefined';
dampMtx = 'undefined';
precompStiffMtx = 'undefined';
precomResVct = 'undefined';

% Title for the output file
title = 'geometrically linear steady-state plane stress analysis';

% Initialize time
t = 0;

% Define tabulation for the output in the command window
tab = '\t';

%% 1. Find the prescribed and the free DOFs of the system

% Prescribed DOFs (DOFs on which either homogeneous or inhomogeneous 
% Dirichlet boundary conditions are prescribed)
prescribedDoFs = mergesorted(homDOFs, inhomDOFs);
prescribedDoFs = unique(prescribedDoFs);

% Free DOFs of the system (actual DOFs over which the solution is computed)
freeDOFs = DOFNumbering;
freeDOFs(ismember(freeDOFs, prescribedDoFs)) = [];

%% 2. Compute the mass matrix of the plate in membrane action
massMtx = computeMassMtxFEMPlateInMembraneAction ...
    (propAnalysis, strMsh, propParameters, propGaussInt);

%% 3. Compute load vector for setting up the initial steady-state computation
if strcmp(outMsg, 'outputEnabled')
    message = 'Computing the load vector for the initial steady-state computation';
    fprintf([tab, '>>', ' ', message]);
end
F = computeLoadVctFEMPlateInMembraneAction ...
    (strMsh, propAnalysis, propNBC, t, propGaussInt, '');
    
%% 4. Solve the initial steady-state problem
if strcmp(outMsg, 'outputEnabled')
    message = 'Performing an initial steady-state computation\n';
    fprintf([tab, '>>', ' ', message]);
end
[dHat, ~, FComplete, minElSize] = solve_FEMNLinearSystem ...
    (propAnalysis, uSaved, uDotSaved, uDDotSaved, strMsh, F, ...
    computeBodyForces, propParameters, dHat, uDot, uDDot, massMtx, dampMtx, ...
    precompStiffMtx, precomResVct, ...
    @computeTangentStiffMtxResVctFEMPlateInMembraneAction, ...
    DOFNumbering, freeDOFs, homDOFs, inhomDOFs, valuesInhomDOFs, ...
    uMeshALE, solve_LinearSystem, propStrDynamics, t, ...
    propNLinearAnalysis, propGaussInt, tab, '');

%% 5. Solve the eigenfrequency problem
if strcmp(outMsg, 'outputEnabled')
    message = '>> Computing the linear stiffness matrix of the structure';
    fprintf([tab, '>>', ' ', message]);
end
[tanMtx, ~, ~] = ...
    computeTangentStiffMtxResVctFEMPlateInMembraneAction ...
    (propAnalysis, dHat, uSaved, uDot, uDotSaved, uMeshALE, ...
    precompStiffMtx, precomResVct, DOFNumbering, strMsh, F, ...
    loadFactor, computeBodyForces, propStrDynamics, t, ...
    propParameters, propGaussInt);

%% 6. Solve the generalized eigenvalue problem to get the eigenmode shapes
if strcmp(outMsg,'outputEnabled')
    message = 'Solving the generalized eigenvalue problem to get the eigenmode shapes\n';
    fprintf([tab, '>>', ' ', message]);
end
[eigenmodeShapesDBC, naturalFrequencies] = ...
    eigs(tanMtx(freeDOFs, freeDOFs), massMtx(freeDOFs, freeDOFs), ...
    numEig, 'sm');
naturalFrequencies = naturalFrequencies*ones(length(naturalFrequencies), 1);
naturalFrequencies = sqrt(naturalFrequencies);
naturalFrequencies = naturalFrequencies./(2*pi);
[noModeShapes, m] = size(eigenmodeShapesDBC);

%% 7. Apply the boundary conditions onto the system
if strcmp(outMsg, 'outputEnabled')
    message = 'Applying the Dirichlet boundary conditions at each mode shape\n\n';
    fprintf([tab, '>>', ' ', message]);
end
eigenmodeShapes = zeros(noModeShapes, m);
for iEig = 1:numEig
    eigenmodeShapes(freeDOFs, iEig) = eigenmodeShapesDBC(:, iEig);
    eigenmodeShapes(homDOFs, iEig) = 0;
    eigenmodeShapes(inhomDOFs, iEig) = 0;
end

%% 8. Write out the results into a file
if isfield(propOutput, 'isOutput')
    if isa(propOutput.isOutput, 'logical')
        if propOutput.isOutput
            if isfield(propOutput, 'writeOutputToFile')
                if isa(propOutput.writeOutputToFile, 'function_handle')
                    fprintf('>> Writting out the results to "%s"\n',strcat(pathToOutput, caseName,'/'));
                    DOF4Output = [1:2:numDOFs - 1
                                  2:2:numDOFs];
                    for iEig = 1:numEig
                        writeOutputFEMPlateInMembraneActionToVTK(propAnalysis, ...
                            propNLinearAnalysis, propStrDynamics, strMsh, ...
                            propParameters, eigenmodeShapes(:, iEig), dDot, ...
                            dDDot, DOF4Output, caseName, pathToOutput, title, ...
                            iEig);
                    end
                else
                    error('Variable propVTK.writeOutputToFile should define a function handle');
                end
            else
                error('Structure propVTK should define variable writeOutputToFile');
            end
        end
    else
        error('Structure propVTK should define boolean isOutput');
    end
end

%% 9. Appendix
if strcmp(outMsg, 'outputEnabled')
    computationalTime = toc;
    fprintf('\nModal analysis took %.2d seconds \n\n', computationalTime);
    fprintf('______________________Nonlinear Analysis Ended____________________\n');
    fprintf('##################################################################\n\n\n');
end

end