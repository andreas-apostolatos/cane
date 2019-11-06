function [dHat,CPHistory,resHistory,hasConverged,minElASize] = ...
    solve_IGAKirchhoffLoveShellNLinear...
    (BSplinePatch,propNLinearAnalysis,solve_LinearSystem,...
    plot_IGANonlinear,graph,outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation 
% 
% Returns the displacement field corresponding to a geometrically
% non-linear anslysis for the isogeometric Kirchhoff-Love shell. 
% Source Reference : 
%
% J. Kiendl "Isogeometric Analysis and Shape Optimal Design for Shell 
% Structures" Ph.D. Thesis, Technische Universtät München (2011)
%
%                Input : 
%         BSplinePatch : Structure containing all the information regarding 
%                        the connections between the multipatches
%  propNLinearAnalysis : Structure on the non-linear analysis settings :
%                                    #Load Steps/Stages#
%                              .noLoadSteps : Selected number of load steps 
%                                             (stages of the loading)
%                                .tolerance : Tolerance for the Newton 
%                                             iterations on the 2-norm
%                                  .maxIter : Maximum number of the Newton 
%                                             iteration
%  solve_LinearSystem : Function handle to the linear equation system
%                       solver
%    plot_IGANonlinear : Function handle to plotting the current
%                        configuation after the solution of each load step
%                graph : On the graphics
%               outMsg : Whether or not to output message on refinement 
%                       progress
%                       'outputEnabled' : enables output information
%
%               Output :
%                 dHat : The displacement field of the patch
%            CPHistory : The deformation history of the Control Points
%           resHistory : The residual history throughout the nonlinear
%                        iterations
%         hasConverged : Flag on whether the nonlinear iterations has
%                        converged or not
%           minElASize : The minimum element area size in the mesh
%
% Function layout :
%
% 0. Read input
%
% 1. Solve the nonlinear system
%
% 2. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('__________________________________________________________________\n');
    fprintf('##################################################################\n');
    fprintf('Static nonlinear analysis for an isogeometric Kirchhoff-Love shell\n');
    fprintf('shell has been initiated\n\n');
    fprintf('Nonlinear scheme : Newton method \n');
    fprintf('Number of load steps = %d \n',propNLinearAnalysis.noLoadSteps);
    fprintf('Residual tolerance = %d \n',propNLinearAnalysis.eps);
    fprintf('Maximum number of nonlinear iterations = %d \n',propNLinearAnalysis.maxIter);
    fprintf('__________________________________________________________________\n\n');

    % start measuring computational time
    tic;
end

%% 0. Read input

% Define analysis type
analysis.type = 'isogemetricKLShellNonlinear';

% Initialize the dummy arrays
dHatSaved = 'undefined';
dHatDotSaved = 'undefined';
dHatDDotSaved = 'undefined';
connections = 'undefined';
dHatDot = 'undefined';
dHatDDot = 'undefined';
computeConstantMatrices = 'undefined';
propCoupling = 'undefined';
massMtx = 'undefined';
dampMtx = 'undefined';
isCosimulationWithEmpire = false;
masterDOFs = [];
slaveDOFs = [];

% Flag on whether the reference configuration is updated
isReferenceUpdated = false;

% Static analysis
propTransientAnalysis.timeDependence = 'steadyState';
propTransientAnalysis.computeProblemMtrcsTransient = 'undefined';
propTransientAnalysis.isStaticStep = true;
t = 0;

% Adjust tabulation
tab = '\t';

% Re-assign the arrays
CP = BSplinePatch.CP;

% Number of Control Points in xi,eta-direction
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));

% Create an element freedom table for the patch in the array
BSplinePatch.DOFNumbering = zeros(nxi,neta,3);
k = 1;
for cpj = 1:neta
    for cpi = 1:nxi
        BSplinePatch.DOFNumbering(cpi,cpj,1) = k;
        BSplinePatch.DOFNumbering(cpi,cpj,2) = k + 1;
        BSplinePatch.DOFNumbering(cpi,cpj,3) = k + 2;

        % Update counter
        k = k + 3;
    end
end

% Create the element freedom table for the BSplinePatch into the array of
% the patches
BSplinePatch.EFTPatches = 1:3*BSplinePatch.noCPs;

% Place the B-Spline patch into an array
BSplinePatches = {BSplinePatch};

% Get number of DOFs
noDOFs = 3*nxi*neta;
BSplinePatches{1}.noDOFs = noDOFs;

% Find the numbering of the DOFs where homogeneous Dirichlet conditions are
% prescribed
homDOFs = BSplinePatch.homDOFs;

% Find the numbering of the free DOFs
freeDOFs = zeros(noDOFs,1);
for i=1:noDOFs
    freeDOFs(i,1) = i;
end
freeDOFs(ismember(freeDOFs,homDOFs)) = [];

% Get the numbering and the values of the DOFs which are prescribed
inhomDOFs = BSplinePatch.inhomDOFs;
valuesInhomDOFs = BSplinePatch.valuesInhomDOFs;

% Initialize the displacement field
dHat = zeros(noDOFs,1);

%% 1. Solve the nonlinear system
[dHat,CPHistory,resHistory,hasConverged,~,~,~,~,BSplinePatches,...
    propCoupling,minElASize] = ...
    solve_IGANLinearSystem...
    (analysis,dHatSaved,dHatDotSaved,dHatDDotSaved,BSplinePatches,connections,dHat,...
    dHatDot,dHatDDot,computeConstantMatrices,massMtx,dampMtx,...
    @computeTangentStiffMtxIGAKirchhoffLoveShellNLinear,...
    @computeUpdatedGeometryIGAThinStructureMultipatches,freeDOFs,homDOFs,...
    inhomDOFs,valuesInhomDOFs,masterDOFs,slaveDOFs,solve_LinearSystem,t,...
    propCoupling,propTransientAnalysis,propNLinearAnalysis,plot_IGANonlinear,...
    isReferenceUpdated,isCosimulationWithEmpire,tab,graph,outMsg);

%% 2. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    fprintf('Static nonlinear analysis took %.2d seconds \n\n',computationalTime);
    fprintf('___________________Static Linear Analysis Ended___________________\n');
    fprintf('##################################################################\n\n\n');
end

end
