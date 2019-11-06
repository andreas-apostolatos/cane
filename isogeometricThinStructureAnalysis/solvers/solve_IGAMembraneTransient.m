function [dHatHistory,resHistory,BSplinePatches,minElAreaSize] = ...
    solve_IGAMembraneTransient...
    (BSplinePatches,computeInitCnds,propNLinearAnalysis,propStrDynamics,...
    propPostproc,solve_LinearSystem,propOutput,pathToOutput,caseName,...
    propEmpireCoSimulation,outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Solves the isogeometric membrane problem and writes out the displacement
% field data into a file to be read by GiD.
%
%                   Input :
%            BSplinePatch : The B-Spline patch array containing:
%                           .p,q : The polynomial orders of the B-Spline 
%                                  surface in both parametric directions
%                        .Xi,Eta : The knot vectors in both parametric 
%                                  directions
%                            .CP : The set of control points and weights
%                       .isNURBS : Flag on whether the basis is a NURBS or 
%                                  a B-Spline
%                    .parameters : Technical parameters for the structure
%                       .homDOFs : The global numbering of the DOFs where 
%                                  homogeneous Dirichlet boundary 
%                                  conditions are applied
%                     .inhomDOFs : The global numbering of the DOFs where
%                                  inhomogeneous Dirichlet boundary 
%                                  conditions are applied
%               .valuesInhomDOFs : The values on the DOFs corresponding to 
%                                  the application of inhomogeneous 
%                                  Dirichlet boundary conditions
%                           .NBC : Structure containing information on the 
%                                  application of the Neumann boundary 
%                                  conditions
%                                       .noCnd : Number of Neumann boundary 
%                                                conditions
%                             .xiLoadExtension : Cell array {.noCnd} 
%                                                containing the load 
%                                                extensions in the xi-
%                                                direction
%                            .etaLoadExtension : Cell array {.noCnd} 
%                                                containing the load 
%                                                extensions in the eta-
%                                                direction
%                               .loadAmplitude : Array (1,.noCnd) 
%                                                containing the load 
%                                                amplitudes
%                               .loadDirection : Array (1,.noCnd) 
%                                                containing the load 
%                                                directions
%                              .computeLoadVct : Cell array {.noCnd} 
%                                                containing the function 
%                                                name for the computation 
%                                                of the load vector
%                              .isConservative : Array (1,.noCnd) of flags 
%                                                indicating whether the 
%                                                load is conservative or 
%                                                not
%         computeInitCnds : Function handle to the computation of the 
%                           initial conditions
%     propNLinearAnalysis : On the properties of the nonlinear solver:
%                           .method : Nonlinear method
%                      .noLoadSteps : Number of load steps
%                              .eps : Solution tolerance on the residual in
%                                     the 2-norm
%                          .maxIter : Maximum number of the nonlinear 
%                                     iterations
%         propStrDynamics : On the properties of the structural dynamics
%                           analysis:
%                   .timeDependence : 'transient' or 'steadyState'
%                           .method : Time integration method
%                           .TStart : Starting time of the simulation
%                             .TEnd : End time of the simulation
%                      .noTimeSteps : Number of time steps
%                               .dt : Time step size
%                          .damping : On the assumed damping:
%                                           .method : Damping method
%                                   .computeDampMtx : Function handle to
%                                                     the computation of 
%                                                     the damping matrix
%                                                     of the system
%     .computeProblemMtrcsTransient : Function handle to the computation of 
%                                     the problem matrices according to the 
%                                     defined time integration method given 
%                                     the system matrix of the steady-state
%                                     system as well as the mass matrix of 
%                                     the problem
%                .computeUpdatedVct : Function handle to the computation of 
%                                     the updated values of the discretized 
%        propPostproc : On the properties of the postprocessing
%                               .resultant : Array of strings containing 
%                                            the name of each resultant to 
%                                            be written out
%                        .computeResultant : Array of strings representing 
%                                            the function handle for the 
%                                            computation of the desirable 
%                                            resultant
%      solve_LinearSystem : Function handle to the linear equation solver
%              propOutput : Properties for writting out the results
%                              .writeOutput : Function handle to outputting 
%                                             the results in a folder
%                           .writeFrequency : Frequency for writting out
%                                             the results
%                            .saveFrequency : Frequency for saving temorary
%                                             results at the time step
%                                             level
%            pathToOutput : Path to the folder where the results are 
%                           written out
%                 caseName : The name of the case in the inputGiD case 
%                            folder
%   propEmpireCoSimulation : Properties for the co-simulation with EMPIRE
%                           .isCoSimulation : Flag on whether co-simulation 
%                                             with EMPIRE is assumed
%                         .isInterfaceLayer : Flag on whether the matlab
%                                             client is used as an
%                                             interface layer
%                              .strMatlabXml : Name of the xml file for 
%                                              connection to Empire
%                   outMsg : Whether or not to output message on refinement 
%                            progress
%                            'outputEnabled' : enables output information
%
%              output :
%         dHatHistory : The history of the Control Point displacements
%          resHistory : The history of the evolution of the residual
%                       throughout the transient analysis
%      BSplinePatches : Updated array of the B-Spline patches
%       minElAreaSize : The minimum element area size in the isogeometric
%                       mesh
%
% Function layout :
%
% 0. Read input
%
% 1. Define function handles
%
% 2. Find the total number of DOFs for the patch including the Lagrange Multipliers and make an EFT for each Lagrange Multipliers field employed
%
% 3. Find the global DOF numbering of the DOFs where homogeneous, inhomogeneous conditions are applied
%
% 4. Solve the transient problem
%
% 5. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('__________________________________________________________________\n');
    fprintf('##################################################################\n');
    fprintf('Transient analysis for an isogeometric membrane has been initiated\n\n');
    if isfield(BSplinePatches{1},'weakDBC')
        if BSplinePatches{1}.weakDBC.noCnd > 0
            fprintf('Weak boundary conditions using the %s method are applied \n',BSplinePatches{1}.weakDBC.method);
            fprintf('over %d boundaries of the B-Spline patch:\n',BSplinePatches{1}.weakDBC.noCnd);
            if strcmp(BSplinePatches{1}.weakDBC.method,'Nitsche')
                if BSplinePatches{1}.weakDBC.estimationStabilPrm == true
                    fprintf('Automatic estimation of the stabilization parameter is enabled \n');
                end
                if isfield(BSplinePatches{1}.weakDBC,'computeConstMtx')
                    if isfield(BSplinePatches{1}.weakDBC,'alpha')
                        fprintf('Manual stabilization parameter chosen as %d\n',BSplinePatches{1}.weakDBC.alpha)
                    else
                        fprintf('Manual stabilization parameter weakDBC.alpha needs to be assigned\n');
                    end
                end
            end
            fprintf('\n');
        end
    end
    if strcmp(propNLinearAnalysis.method,'undefined')
        fprintf('Geometrically linear analysis chosen \n');
    else
        fprintf('Geometrically nonlinear analysis using the %s method \n',propNLinearAnalysis.method);
        fprintf('chosen:\n');
    end
    fprintf('Number of load steps = %d \n',propNLinearAnalysis.noLoadSteps);
    fprintf('Residual tolerance = %d \n',propNLinearAnalysis.eps);
    fprintf('Maximum number of nonlinear iterations = %d \n\n',propNLinearAnalysis.maxIter);
    if propEmpireCoSimulation.isCoSimulation
        fprintf('Co-simulation through Empire \n');
        fprintf('---------------------------- \n\n');
        if isempty(propEmpireCoSimulation.strMatlabXml)
            error('No xml file for the Matlab code is provided');
        else
            fprintf('Matlab xml file for connection to Empire: %s.xml\n',propEmpireCoSimulation.strMatlabXml);
        end
        fprintf('\n');
    end
    fprintf('__________________________________________________________________\n\n');

    % start measuring computational time
    tic;
end

%% 0. Read input

% Define the analysis type
analysis.type = 'isogeometricMembraneAnalysis';

% Title for the output file
if ~strcmp(propNLinearAnalysis.method,'undefined')
    title = 'Geometrically nonlinear transient isogeometric membrane analysis';
else
    title = 'Geometrically linear transient isogeometric membrane analysis';
end

% Initialize the dummy arrays
propCoupling = 'undefined';
nodesALE = 'undefined';
connections = 'undefined';
computeUpdatedMesh = 'undefined';

% Flag on whether the reference configuration is updated
isReferenceUpdated = false;

% On the application of weak Dirichlet boundary conditions
isWeakDBC = false;
if ~isempty(BSplinePatches{1}.weakDBC)
    if isfield(BSplinePatches{1}.weakDBC,'noCnd')
        if BSplinePatches{1}.weakDBC.noCnd > 0
            isWeakDBC = true;
        end
    end
end

% Adjust tabulation
tab = '\t';

%% 1. Define function handles

% Function handle to the computation of the tangent stiffness matrix and
% residual vector for a single patch geometry
if BSplinePatches{1}.noElmnts ~= 1
    computeProblemMtrcsSteadyState = @computeTangentStiffMtxResVctIGAMembraneNLinear;
else
    computeProblemMtrcsSteadyState = @computeTangentStiffMtxResVctIGAMembraneNLinearOutdated;
end

% Function handle to the computation of the mass matrix
computeMassMtx = @computeIGAMassMtxThinStructure;

% Function handle to the computation of the updated deformed geometry 
computeUpdatedGeometry = @computeUpdatedGeometryIGAThinStructureMultipatches;

% Compute constant problem matrices in case of the application of weak
% boundary conditions
if isWeakDBC
    if isfield(BSplinePatches{1}.weakDBC,'method')
        if strcmp(BSplinePatches{1}.weakDBC.method,'penalty') || ...
                (strcmp(BSplinePatches{1}.weakDBC.method,'nitsche') && ...
                ~BSplinePatches{1}.weakDBC.estimationStabilPrm)
            computeConstantProblemMatrices = @computeWeakDBCMtxPenaltyIGAMembrane;
        elseif strcmp(BSplinePatches{1}.weakDBC.method,'lagrangeMultipliers')
            computeConstantProblemMatrices = @computeWeakDBCMtxLagrangeMultipliersIGAMembrane;
        elseif strcmp(BSplinePatches{1}.weakDBC.method,'nitsche') && ...
                BSplinePatches{1}.weakDBC.estimationStabilPrm
            computeConstantProblemMatrices = 'undefined';
        elseif ~strcmp(BSplinePatches{1}.weakDBC.method,'penalty') && ...
                ~strcmp(BSplinePatches{1}.weakDBC.method,'lagrangeMultipliers') && ...
                ~strcmp(BSplinePatches{1}.weakDBC.method,'nitsche')
            error('Define a valid method in BSplinePatches{1}.weakDBC.method');
        end
    else
        error('Define a method in BSplinePatches{1}.weakDBC');
    end
else
    computeConstantProblemMatrices = 'undefined';
end

%% 2. Find the total number of DOFs for the patch including the Lagrange Multipliers and make an EFT for each Lagrange Multipliers field employed

% Assign an EFT for each Lagrange Multipliers field
noDOFsPatchLM = 0;
if isWeakDBC
    if strcmp(BSplinePatches{1}.weakDBC.method,'lagrangeMultipliers')
        for iCnd = 1:BSplinePatches{1}.weakDBC.noCnd
            noDOFsPatchLMCnd = 3*length(BSplinePatches{1}.weakDBC.lambda{iCnd}.CP(:,1));
            noDOFsPatchLM = noDOFsPatchLM + noDOFsPatchLMCnd;
            if iCnd == 1
                index = 3*BSplinePatches{1}.noCPs;
            else
                index = BSplinePatches{1}.weakDBC.lambda{iCnd-1}.EFT(length(BSplinePatches{1}.weakDBC.lambda{iCnd-1}.EFT));
            end
            BSplinePatches{1}.weakDBC.lambda{iCnd}.EFT = index + 1:index + noDOFsPatchLMCnd;
        end
    end
end
BSplinePatches{1}.noDOFs = 3*BSplinePatches{1}.noCPs + noDOFsPatchLM;
if propEmpireCoSimulation.isCoSimulation
    BSplinePatches{1}.noDOFsEmpire = 3*BSplinePatches{1}.noCPs;
end

% Compute the number of DOFs
noDOFs = BSplinePatches{1}.noDOFs;

% Create the element freedom table for the BSplinePatch into the array of
% the patches
BSplinePatches{1}.EFTPatches = 1:BSplinePatches{1}.noDOFs;

% Create a DOF numbering for each Lagrange Multipliers field of the patch
if isWeakDBC
    if strcmp(BSplinePatches{1}.weakDBC.method,'lagrangeMultipliers')
        for iCnd = 1:BSplinePatches{1}.weakDBC.noCnd
            % Get the number of Control Points
            nxiLambda = length(BSplinePatches{1}.weakDBC.lambda{iCnd}.Xi);

            % Initialize the field of the DOF numbering
            BSplinePatches{1}.weakDBC.lambda{iCnd}.DOFNumbering = zeros(nxiLambda,3);

             % Compute the entries of the DOF numbering array
            k = 1;
            for cpi = 1:nxiLambda
                BSplinePatches{1}.weakDBC.lambda{iCnd}.DOFNumbering(cpi,1) = k;
                BSplinePatches{1}.weakDBC.lambda{iCnd}.DOFNumbering(cpi,2) = k + 1;
                BSplinePatches{1}.weakDBC.lambda{iCnd}.DOFNumbering(cpi,3) = k + 2;

                % Update counter
                k = k + 3;
            end
        end
    end
end

%% 3. Find the global DOF numbering of the DOFs where homogeneous, inhomogeneous conditions are applied

% Get the numbering and the values of the DOFs which are prescribed
inhomDOFs = BSplinePatches{1}.inhomDOFs;
valuesInhomDOFs = BSplinePatches{1}.valuesInhomDOFs;

% Find the numbering of the DOFs where homogeneous Dirichlet conditions are
% prescribed
homDOFs = BSplinePatches{1}.homDOFs;

% Clear the inhomogeneous Dirichlet boundary conditions from the
% homogeneous
homDOFs(ismember(homDOFs,inhomDOFs)) = [];

% Find the numbering of the free DOFs
freeDOFs = zeros(noDOFs,1);
for i = 1:noDOFs
    freeDOFs(i,1) = i;
end
freeDOFs(ismember(freeDOFs,homDOFs)) = [];
freeDOFs(ismember(freeDOFs,inhomDOFs)) = [];

% Master and slave DOFs
masterDOFs = [];
slaveDOFs = [];

%% 4. Solve the transient problem
[dHatHistory,resHistory,BSplinePatches,~,minElAreaSize] = ...
    solve_IGATransientAnalysis...
    (analysis,BSplinePatches,connections,freeDOFs,homDOFs,inhomDOFs,...
    valuesInhomDOFs,masterDOFs,slaveDOFs,nodesALE,computeInitCnds,...
    @solve_IGANLinearSystem,computeConstantProblemMatrices,...
    computeMassMtx,computeProblemMtrcsSteadyState,computeUpdatedMesh,...
    computeUpdatedGeometry,solve_LinearSystem,propCoupling,...
    propStrDynamics,propNLinearAnalysis,propPostproc,caseName,...
    pathToOutput,title,propOutput,isReferenceUpdated,...
    propEmpireCoSimulation,tab,outMsg);

%% 5. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    fprintf('Transient nonlinear analysis took %.2d seconds \n\n',computationalTime);
    fprintf('________________Transient Nonlinear Analysis Ended________________\n');
    fprintf('##################################################################\n\n\n');
end

end

