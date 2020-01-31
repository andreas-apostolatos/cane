function [x_Min] = computeStructureBoundary(msh,propALE)

%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Matthew Keller
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
% 1. Save the coordinates of the nodes before the ALE motion
%
% 2. Loop over all the ALE nodes and assign the homogeneous and inhomogeneous Dirichlet boundary conditions for the ALE nodes 
% ->
%    2i. Find the node ID in the node cloud
%
%   2ii. Get the coordinates of the node
%
% 3. Return the minimum x value

%% 0 Initialize structure array counter
i = 1;

%% The following steps are only performed if there exists an ALE boundary
if ~isempty(propALE)
    
    %% 2. Loop over all the ALE nodes and assign coordinates based on function handle designation
    for counterALE = 1:length(propALE.nodes(:,1))
        
        % Exacute if the node function handle equates to the structure handle
        if strcmp(propALE.fctHandle((counterALE),:),'computeALEMM')
            %% 2i. Find the node ID in the node cloud
            nodeID = propALE.nodes(counterALE,1);

            %% 2ii. Get the coordinates of the node
            x(i) = msh.initialNodes(nodeID,1);
            y(i) = msh.initialNodes(nodeID,2);
            z(i) = msh.initialNodes(nodeID,3);
            
            % Iterate counter
            i = i+1;
        end
    end
    
    %% 3. Return the minimum x coordinate
    x_Min = min(x);
    
end