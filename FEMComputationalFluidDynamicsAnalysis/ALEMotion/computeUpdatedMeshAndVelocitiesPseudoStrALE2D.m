function [msh, uMeshALE, inhomDOFs, valuesInhomDOFs] = ...
    computeUpdatedMeshAndVelocitiesPseudoStrALE2D ...
    (msh, homDOFs, inhomDOFs, valuesInhomDOFs, nodesSaved, propALE, ...
    solve_LinearSystem, propTransientAnalysis, t)
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
%             inhomDOFs : The updated vector with the global numbering of 
%                         the DOFs where inhomogeneous boundary conditions 
%                         are applied
%       valuesInhomDOFs : the updated vector with the values of the 
%                         inhomogeneous Dirichlet boundary conditions at 
%                         each node
%
% Function layout :
%
% 0. Read input
%
% The following steps are only performed if there exists an ALE boundary
%
% 1. Loop over all the ALE nodes and assign the homogeneous and inhomogeneous Dirichlet boundary conditions for the ALE nodes 
% ->
%    1i. Find the node ID in the node cloud
%
%   1ii. Get the coordinates of the node
%
%  1iii. Get function handle for ALE motion computation
%
%   1iv. Compute the motion at the given node
%
%    1v. Assign the homogeneous or inhomogeneous boundary condition for the given node
% <-
%
% 2. Get the array of the ALE nodes with a prescribed motion
%
% 3. Get the free DOFs of the system
%
% 4. Solve the pseudostructural problem for the ALE motion
%
% 5. Move the nodes on the mesh following the prescribed motion
%
% 6. Loop over all the ALE nodes, compute the mesh velocity and updated the inhomogeneous Dirichlet boundary conditions
% ->
%    6i. Find the node ID
%
%   6ii. Compute the velocity of the node using a first order interpolation
%
%  6iii. Collect all the inhomogeneous boundary conditions from the ALE motion into an array
% <-
%
% 7. Collect all the DOFs where inhomogeneous boundary conditions are prescribed into an array
%
% 8. Compute the mesh velocity at the interior nodes of mesh using a first order interpolation
%
%% Function main body

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
computeBodyForces = @(x,y,z,t) zeros(length(x),1,3);

% Initialize boundary conditions for the ALE motion
homDBCALE = [];
inhomDBCALE = [];
valuesInhomDBCALE = [];
inhomDOFsALE = [];
valuesInhomDOFsALE = [];
homDOFsALE = [];

% Number of nodes in the domain
noNodes = length(msh.nodes(:,1));

% Number of DOFs in the domain
noDOFs = 2*noNodes;

% Initialize the array of the velocity of the nodes in the mesh due to the
% ALE motion
uMeshALE = zeros(3*noNodes,1);

% Make a global DOF numbering for the ALE problem
DOFNumberingALE = 1:noDOFs;

% Initialize vector of DOFs for the ALE problem
dHatALE = zeros(noDOFs,1);

% Initialize counters
counterIhnHomDOFs = 1;
counterHomDBC = 1;
counterInhomDBC = 1;

%% The following steps are only performed if there exists an ALE boundary
if ~isempty(propALE)
    %% 1. Loop over all the ALE nodes and assign the homogeneous and inhomogeneous Dirichlet boundary conditions for the ALE nodes
    for iALE = 1:length(propALE.nodes(:,1))
        %% 1i. Find the node ID in the node cloud
        nodeID = propALE.nodes(iALE,1);

        %% 1ii. Get the coordinates of the node
        x = msh.initialNodes(nodeID,1);
        y = msh.initialNodes(nodeID,2);
        z = msh.initialNodes(nodeID,3);

        %% 1iii. Get function handle for ALE motion computation
        computeALEMotion = str2func(propALE.fctHandle((iALE),:));

        %% 1iv. Compute the motion at the given node
        [dx,dy,~] = computeALEMotion(x, y, z, t, propUser);
        
        %% 1v. Assign the homogeneous or inhomogeneous boundary condition for the given node
        if dx == 0
            homDBCALE(1,counterHomDBC) = 2*nodeID-1;
            counterHomDBC = counterHomDBC + 1;
        elseif isnan(dx)
            % Nothing is done if dx is not a number
        else
            inhomDBCALE(1,counterInhomDBC) = 2*nodeID-1;
            valuesInhomDBCALE(1,counterInhomDBC) = dx;
            counterInhomDBC = counterInhomDBC + 1;
        end
        if dy == 0
            homDBCALE(1,counterHomDBC) = 2*nodeID;
            counterHomDBC = counterHomDBC + 1;
        elseif isnan(dy)
            % Nothing is done if dy is not a number
        else
            inhomDBCALE(1,counterInhomDBC) = 2*nodeID;
            valuesInhomDBCALE(1,counterInhomDBC) = dy;
            counterInhomDBC = counterInhomDBC + 1;
        end
    end

    %% 2. Get the array of the ALE nodes with a prescribed motion
    prescribedDOFsALE = mergesorted(homDBCALE,inhomDBCALE);
    prescribedDOFsALE = unique(prescribedDOFsALE);

    %% 3. Get the free DOFs of the system
    freeDOFsALE = DOFNumberingALE;
    freeDOFsALE(prescribedDOFsALE) = [];

     %% 4. Solve the pseudostructural problem for the ALE motion
    [dHatALE,~,~,~] = solve_FEMLinearSystem...
        (propAnalysisALE, uSavedALE, uDotSavedALE, uDDotSavedALE, msh,...
        zeros(2*length(msh.nodes(:,1)),1), computeBodyForces, ...
        parametersALE, dHatALE, uDotALE, uDDotALE, massMtx, dampMtx, ...
        precompStiffMtx, precomResVct, ...
        @computeStiffMtxAndLoadVctFEMPlateInMembraneActionCST, ...
        DOFNumberingALE, freeDOFsALE, homDBCALE, inhomDBCALE, ...
        valuesInhomDBCALE, uMeshALE, solve_LinearSystem, propStrDynamics, ...
        propNLinearAnalysis, propIntALE, tab, '');

    %% 5. Move the nodes on the mesh following the prescribed motion
    for i = 1:noNodes
        msh.nodes(i,1) = msh.initialNodes(i,1) + dHatALE(2*i-1);
        msh.nodes(i,2) = msh.initialNodes(i,2) + dHatALE(2*i);
    end
    
    %% Debugging
%     graph.index = 1;
%     graph.index = plot_referenceConfigurationFEMPlateInMembraneAction...
%         (msh,'undefined',zeros(length(msh.nodes(:,1)),1),[],graph,'');
%     graph.index = graph.index + 1;

    %% 6. Loop over all the ALE nodes, compute the mesh velocity and updated the inhomogeneous Dirichlet boundary conditions
    for i = 1:length(propALE.nodes(:,1))
        %% 6i. Find the node ID
        nodeID = propALE.nodes(i,1);
        
        %% 6ii. Compute the velocity of the node using a first order interpolation
        if strcmp(propTransientAnalysis.timeDependence,'TRANSIENT')
            ux = (msh.nodes(nodeID,1) - nodesSaved(nodeID,1))/propTransientAnalysis.dt;
            uy = (msh.nodes(nodeID,2) - nodesSaved(nodeID,2))/propTransientAnalysis.dt;
        elseif strcmp(propTransientAnalysis.timeDependence,'STEADY_STATE')
            ux = 0;
            uy = 0;
        else
            error('Neither steady-state nor transient analysis is selected, see input file');
        end
        
        %% 6iii. Collect all the inhomogeneous boundary conditions from the ALE motion into an array
        
        % Velocity in x-direction
        if  ux ~= 0 && strcmp(propALE.fctHandle(i,:), 'computeALEMF')
            inhomDOFsALE(1,counterInhomDOFs) = 3*nodeID - 2;
            valuesInhomDOFsALE(1,counterInhomDOFs) = ux;
            counterInhomDOFs = counterInhomDOFs + 1;
        end
        
        % Velocity in y-direction
        if uy ~= 0 && strcmp(propALE.fctHandle(i,:), 'computeALEMF')
            inhomDOFsALE(1,counterIhnHomDOFs) = 3*nodeID - 1;
            valuesInhomDOFsALE(1,counterIhnHomDOFs) = uy;
            counterIhnHomDOFs = counterIhnHomDOFs + 1;
        end
    end
    
    %% 7. Collect all the DOFs where inhomogeneous boundary conditions are prescribed into an array
    inhomDOFs = horzcat(inhomDOFs, inhomDOFsALE);
    valuesInhomDOFs = horzcat(valuesInhomDOFs, valuesInhomDOFsALE);
    [inhomDOFs, indexSorting] = sort(inhomDOFs);
    valuesInhomDOFs = valuesInhomDOFs(indexSorting);

    %% 8. Compute the mesh velocity at the interior nodes of mesh using a first order interpolation
    for i = 1:length(msh.nodes(:, 1))
        if strcmp(propTransientAnalysis.timeDependence, 'TRANSIENT')
            ux = (msh.nodes(i, 1) - nodesSaved(i, 1))/propTransientAnalysis.dt;
            uy = (msh.nodes(i, 2) - nodesSaved(i, 2))/propTransientAnalysis.dt;
        elseif strcmp(propTransientAnalysis.timeDependence, 'STEADY_STATE')
            ux = 0;
            uy = 0;
        else
            error('Neither steady-state nor transient analysis is selected, see input file');
        end
        uMeshALE(3*i - 2, 1) = ux;
        uMeshALE(3*i - 1, 1) = uy;
        uMeshALE(3*i, 1) = 0;
    end
end

end
