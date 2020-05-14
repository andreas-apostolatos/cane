function [msh, uMeshALE, homDOFs, inhomDOFs, valuesInhomDOFs, freeDOFs] = ...
    computeUpdatedMeshAndVelocitiesPseudoStrALE2D ...
    (msh, homDOFs, inhomDOFs, valuesInhomDOFs, freeDOFs, nodesSaved, ...
    propALE, solve_LinearSystem, propTransientAnalysis, t)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the updated mesh node coordinates and the according velocities on
% the moving boundaries. For the mesh motion the linear pseudo-structural
% solver is employed.
%
%                 Input :
%                   msh : Nodes and elements in the low order mesh:
%                               .nodes : The nodes in the mesh
%                            .elements : The elements in the mesh
%               homDOFs : The global numbering of the DOFs where 
%                         homogeneous boundary conditions are applied
%             inhomDOFs : The global numbering of the DOFs where 
%                         inhomogeneous boundary conditions are applied
%       valuesInhomDOFs : The values of the inhomogeneous Dirichlet 
%                         boundary conditions at each node
%              freeDOFs : The global numbering of the free DOFs to be 
%                         solved for
%            nodesSaved : Saved nodal coordinates from the previous time
%                         step
%               propALE : Properties regarding the ALE boundary
%                           .nodes : The sequence of the nodal coordinates
%                                    on the ALE boundary
%                       .fcthandle : Function handle to the computation of
%                                    the ALE motion
%                         propUser : Extra user-defined parameters
%    solve_LinearSystem : Function handle to the solver for the linear 
%                         equation system
% propTransientAnalysis : On the transient analysis :
%                               .method : The time integration method
%                            .alphaBeta : (parameter for the Bossak scheme)
%                                .gamma : (parameter for the Bossak scheme)
%                               .TStart : Start time of the simulation
%                                 .TEnd : End time of the simulation
%                                   .nT : Number of time steps
%                                   .dt : Time step
%                     t : The current time instance
%
%                Output :
%                   msh : The mesh with the updated nodal coordinates
%              uMeshALE : The velocity of the nodes in the moving mesh due 
%                         to the moving boundary
%               homDOFs : The updated vector with the global numbering of 
%                         the DOFs where homogeneous boundary conditions 
%                         are applied
%             inhomDOFs : The updated vector with the global numbering of 
%                         the DOFs where inhomogeneous boundary conditions 
%                         are applied
%       valuesInhomDOFs : the updated vector with the values of the 
%                         inhomogeneous Dirichlet boundary conditions at 
%                         each node
%              freeDOFs : The updated vector of the free DOFs to be solved 
%                         for
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the ALE nodes to obtain the boundary conditions and set up the mesh motion problem
% ->
%    1i. Find the node ID in the node cloud
%
%   1ii. Get the coordinates of the node
%
%  1iii. Get function handle for ALE motion computation
%
%   1iv. Compute the motion at the given node
% <-
%
% 2. Get the array of the ALE nodes with a prescribed motion
%
% 3. Get the free DOFs of the system for the ALE problem
%
% 4. Solve the pseudostructural problem for the ALE motion
%
% 5. Move the nodes on the mesh following the prescribed motion
%
% 6. Loop over all the ALE nodes and adjust the boundary conditions for the fluid problem
% ->
%    6i. Find the node ID
%
%   6ii. Check whether the node belongs to a node on the free stream
%
%  6iii. Compute the velocity of the node using a first order interpolation
%
%   6iv. Collect all the inhomogeneous boundary conditions from the ALE motion into an array
% <-
%
% 7. Collect all the DOFs where homogeneous boundary conditions are prescribed into an array
%
% 8. Collect all the DOFs where inhomogeneous boundary conditions are prescribed into an array
%
% 9. Get the array of the fluid nodes with a prescribed motion
%
% 10. Get the free DOFs of the system for the fluid problem
%
% 11. Compute the mesh velocity at the interior nodes of mesh using a first order interpolation
%
%% Function main body
if isempty(propALE) || ischar(propALE)
    error('propALE is undefined and mesh motion is called');
else
    if ~isfield(propALE, 'nodes')
        error('propALE must define a field with name nodes')
    end
    if ~isfield(propALE, 'fctHandle')
        error('propALE must define a field with name fctHandle')
    end
    if ~isfield(propALE, 'isFree')
        error('propALE must define a field with name isFree')
    end
    if length(propALE.nodes(:, 1)) ~= length(propALE.fctHandle(:, 1)) || ...
            length(propALE.nodes(:, 1)) ~= length(propALE.isFree(:, 1)) || ...
            length(propALE.fctHandle(:, 1)) ~= length(propALE.isFree(:, 1))
        error('Arrays propALE.nodes, propALE.fctHandle and propALE.isFree must have the same length');
    end
end

%% 0. Read input

% Initialize dummy arrays
uSavedALE = 'undefined'; 
uDotSavedALE = 'undefined'; 
uDDotSavedALE = 'undefined'; 
uDotALE = 'undefined';
uDDotALE = 'undefined';
massMtx = 'undefined';
dampMtx = 'undefined';
precompStiffMtx = 'undefined';
precomResVct = 'undefined';
propNLinearAnalysis = 'undefined';
propStrDynamics = 'undefined';
tab = 'undefined';
propIntALE.type = 'default';
propAnalysisALE.type = 'planeStress';

% The material properties of the pseudo-structural solver
parametersALE.nue = 0;
parametersALE.E = 1e3;

% User-defined properties
if isfield(propALE,'propUser')
    propUser = propALE.propUser;
else
    propUser = 'undefined';
end

% Zero body forces
computeBodyForces = @(x,y,z,t) zeros(length(x), 1, 3);

% Initialize vectors containing the global numbering of DOFs which are 
% either constrained or free regarding the pseudo-structural solver
homDOFsALE_pseudoStr = [];
inhomDOFsALE_preudoStr = [];
valuesInhomDOFsALE_pseudoStr = [];

% Initialize vectors containing the global numbering of DOFs which are 
% either constrained or free regarding the fluid solver after moving the
% mesh
homDOFsALE = [];
inhomDOFsALE = [];
valuesInhomDOFsALE = [];

% Number of nodes in the domain
numNodes = length(msh.nodes(:,1));

% Number of DOFs in the domain
numDOFs = 2*numNodes;

% Initialize the array of the velocity of the nodes in the mesh due to the
% ALE motion
uMeshALE = zeros(3*numNodes,1);

% Make a global DOF numbering for the ALE problem
DOFNumberingALE_presudoStr = 1:numDOFs;

% Initialize vector of DOFs for the ALE problem
dHatALE = zeros(numDOFs, 1);

%% 1. Loop over all the ALE nodes to obtain the boundary conditions and set up the mesh motion problem
for iALE = 1:length(propALE.nodes(:,1))
    %% 1i. Find the node ID in the node cloud
    nodeID = propALE.nodes(iALE,1);

    %% 1ii. Get the coordinates of the node
    x = msh.initialNodes(nodeID,2);
    y = msh.initialNodes(nodeID,3);
    z = msh.initialNodes(nodeID,4);

    %% 1iii. Get function handle for ALE motion computation
    computeALEMotion = str2func(propALE.fctHandle{iALE, 1});

    %% 1iv. Compute the motion at the given node
    [dx, dy, ~] = computeALEMotion(x, y, z, t, propUser);

    %% 1v. Assign the homogeneous or inhomogeneous boundary condition for the given node
    if dx == 0
        homDOFsALE_pseudoStr = horzcat(homDOFsALE_pseudoStr, 2*nodeID - 1);
    else
        if ~isnan(dx)
            inhomDOFsALE_preudoStr = horzcat(inhomDOFsALE_preudoStr, 2*nodeID - 1);
            valuesInhomDOFsALE_pseudoStr = horzcat(valuesInhomDOFsALE_pseudoStr, dx);
        end
    end
    if dy == 0
        homDOFsALE_pseudoStr = horzcat(homDOFsALE_pseudoStr, 2*nodeID);
    else
        if ~isnan(dy)
            inhomDOFsALE_preudoStr = horzcat(inhomDOFsALE_preudoStr, 2*nodeID);
            valuesInhomDOFsALE_pseudoStr = horzcat(valuesInhomDOFsALE_pseudoStr, dy);
        end
    end
end

%% 2. Get the array of the ALE nodes with a prescribed motion
prescribedDOFsALE_pseudoStr = mergesorted(homDOFsALE_pseudoStr, inhomDOFsALE_preudoStr);
prescribedDOFsALE_pseudoStr = unique(prescribedDOFsALE_pseudoStr);

%% 3. Get the free DOFs of the system for the ALE problem
freeDOFsALE_pseudoStr = DOFNumberingALE_presudoStr;
freeDOFsALE_pseudoStr(ismember(freeDOFsALE_pseudoStr, prescribedDOFsALE_pseudoStr)) = [];

 %% 4. Solve the pseudostructural problem for the ALE motion
[dHatALE,~,~,~] = solve_FEMLinearSystem...
    (propAnalysisALE, uSavedALE, uDotSavedALE, uDDotSavedALE, msh,...
    zeros(2*length(msh.nodes(:, 1)), 1), computeBodyForces, ...
    parametersALE, dHatALE, uDotALE, uDDotALE, massMtx, dampMtx, ...
    precompStiffMtx, precomResVct, ...
    @computeStiffMtxAndLoadVctFEMPlateInMembraneActionCST, ...
    DOFNumberingALE_presudoStr, freeDOFsALE_pseudoStr, homDOFsALE_pseudoStr, inhomDOFsALE_preudoStr, ...
    valuesInhomDOFsALE_pseudoStr, uMeshALE, solve_LinearSystem, propStrDynamics, ...
    t, propNLinearAnalysis, propIntALE, tab, '');

%% 5. Move the nodes on the mesh following the prescribed motion
msh.nodes(:, 2) = msh.initialNodes(:, 2) + dHatALE(1:2:end - 1, 1);
msh.nodes(:, 3) = msh.initialNodes(:, 3) + dHatALE(2:2:end, 1);

%% Debugging
%     graph.index = 1;
%     graph.index = plot_referenceConfigurationFEMPlateInMembraneAction...
%         (msh,'undefined',zeros(length(msh.nodes(:,1)),1),[],graph,'');
%     graph.index = graph.index + 1;

%% 6. Loop over all the ALE nodes and adjust the boundary conditions for the fluid problem
for iALE = 1:length(propALE.nodes(:, 1))
    %% 6i. Find the node ID
    nodeID = propALE.nodes(iALE, 1);
    
    %% 6ii. Check whether the node belongs to a node on the free stream
    isFree = propALE.isFree(iALE, 1);

    %% 6iii. Compute the velocity of the node using a first order interpolation
    if strcmp(propTransientAnalysis.timeDependence, 'TRANSIENT')
        ux = (msh.nodes(nodeID,2) - nodesSaved(nodeID,2))/propTransientAnalysis.dt;
        uy = (msh.nodes(nodeID,3) - nodesSaved(nodeID,3))/propTransientAnalysis.dt;
    elseif strcmp(propTransientAnalysis.timeDependence, 'STEADY_STATE')
        ux = 0;
        uy = 0;
    else
        error('Neither steady-state nor transient analysis is selected, see input file');
    end

    %% 6iv. Collect all the inhomogeneous boundary conditions from the ALE motion into an array
    if ~isFree
        if ux ~= 0
            inhomDOFsALE = horzcat(inhomDOFsALE, 3*nodeID - 2);
            valuesInhomDOFsALE = horzcat(valuesInhomDOFsALE, ux);
        else
            homDOFsALE = horzcat(homDOFsALE, 3*nodeID - 2);
        end
        if uy ~= 0
            inhomDOFsALE = horzcat(inhomDOFsALE, 3*nodeID - 1);
            valuesInhomDOFsALE = horzcat(valuesInhomDOFsALE, uy);
        else
            homDOFsALE = horzcat(homDOFsALE, 3*nodeID - 1);
        end
    end
end

%% 7. Collect all the DOFs where homogeneous boundary conditions are prescribed into an array
homDOFs = horzcat(homDOFs, homDOFsALE);
[homDOFs, ~] = sort(homDOFs);

%% 8. Collect all the DOFs where inhomogeneous boundary conditions are prescribed into an array
inhomDOFs = horzcat(inhomDOFs, inhomDOFsALE);
valuesInhomDOFs = horzcat(valuesInhomDOFs, valuesInhomDOFsALE);
[inhomDOFs, indexSorting] = sort(inhomDOFs);
valuesInhomDOFs = valuesInhomDOFs(indexSorting);

%% 9. Get the array of the fluid nodes with a prescribed motion
prescribedDOFs = mergesorted(homDOFs, inhomDOFs);
prescribedDOFs = unique(prescribedDOFs);

%% 10. Get the free DOFs of the system for the fluid problem
freeDOFs(ismember(freeDOFs, prescribedDOFs)) = [];

%% 11. Compute the mesh velocity at the interior nodes of mesh using a first order interpolation
for iALE = 1:length(msh.nodes(:, 1))
    if strcmp(propTransientAnalysis.timeDependence, 'TRANSIENT')
        ux = (msh.nodes(iALE, 2) - nodesSaved(iALE, 2))/propTransientAnalysis.dt;
        uy = (msh.nodes(iALE, 3) - nodesSaved(iALE, 3))/propTransientAnalysis.dt;
    elseif strcmp(propTransientAnalysis.timeDependence, 'STEADY_STATE')
        ux = 0;
        uy = 0;
    else
        error('Neither steady-state nor transient analysis is selected, see input file');
    end
    uMeshALE(3*iALE - 2, 1) = ux;
    uMeshALE(3*iALE - 1, 1) = uy;
    uMeshALE(3*iALE, 1) = 0;
end

end
