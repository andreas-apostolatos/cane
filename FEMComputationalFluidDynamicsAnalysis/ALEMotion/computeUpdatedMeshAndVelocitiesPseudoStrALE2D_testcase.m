function [msh,uMeshALE,inhomDOFs,valuesInhomDOFs] = computeUpdatedMeshAndVelocitiesPseudoStrALE2D_testcase...
    (msh,homDOFs,inhomDOFs,valuesInhomDOFs,nodesALE,solve_LinearSystem,propTransientAnalysis,t)
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
%             inhomDOFs : The global numbering of the DOFs where 
%                         inhomogeneous boundary conditions are applied
%       valuesInhomDOFs : The values of the inhomogeneous Dirichlet 
%                         boundary conditions at each node
%              nodesALE : The nodes on the ALE boundary:
%                           .nodes : The sequence of the nodal coordinates
%                                    on the ALE boundary
%                       .fcthandle : Function handle to the computation of
%                                    the ALE motion
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
% 1. Save the coordinates of the nodes before the ALE motion
%
% 2. Loop over all the ALE nodes and assign the homogeneous and inhomogeneous Dirichlet boundary conditions for the ALE nodes 
% ->
%    2i. Find the node ID in the node cloud
%
%   2ii. Get the coordinates of the node
%
%  2iii. Get function handle for ALE motion computation
%
%   2iv. Compute the motion at the given node
%
%    2v. Assign the homogeneous or inhomogeneous boundary condition for the given node
% <-
%
% 3. Get the array of the ALE nodes with a prescribed motion
%
% 4. Get the free DOFs of the system
%
% 5. Solve the pseudostructural problem for the ALE motion
%
% 6. Move the nodes on the mesh following the prescribed motion
%
% 7. Loop over all the ALE nodes, compute the mesh velocity and updated the inhomogeneous Dirichlet boundary conditions
% ->
%    7i. Find the node ID
%
%   7ii. Compute the velocity of the node using a first order interpolation
%
%  7iii. Collect all the inhomogeneous boundary conditions from the ALE motion into an array
% <-
%
% 8. Collect all the DOFs where inhomogeneous boundary conditions are prescribed into an array
%
% 9. Compute the mesh velocity at the interior nodes of mesh using a first order interpolation
%
%% Function main body

%% 0. Read input

% Number of nodes in the domain
noNodes = length(msh.nodes(:,1));

% Initialize the array of the velocity of the nodes in the mesh due to the
% ALE motion
uMeshALE = zeros(3*noNodes,1);

%% The following steps are only performed if there exists an ALE boundary
if ~isempty(nodesALE)
    %% 1. Save the coordinates of the nodes before the ALE motion
    nodesSaved = msh.nodes;
    
    %% 2. Loop over all the ALE nodes and assign the homogeneous and inhomogeneous Dirichlet boundary conditions for the ALE nodes
    for counterALE = 1:length(nodesALE.nodes(:,1))
        %% 2i. Find the node ID in the node cloud
        nodeID = nodesALE.nodes(counterALE,1);

        %% 2ii. Get the coordinates of the node
        x = msh.initialNodes(nodeID,1);
        y = msh.initialNodes(nodeID,2);
        z = msh.initialNodes(nodeID,3);

        %% 2iii. Get function handle for ALE motion computation
        computeALEMotion = str2func(nodesALE.fctHandle((counterALE),:));

        %% 2iv. Compute the motion at the given node
        [dx,dy,~] = computeALEMotion(x,y,z,t);
        
        %% 2v. Assign the homogeneous or inhomogeneous boundary condition for the given node
        msh.nodes(nodeID,1) = msh.initialNodes(nodeID,1) + dx;
        msh.nodes(nodeID,2) = msh.initialNodes(nodeID,2) + dy;
    end
end