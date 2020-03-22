function [x_Base_Min,x_Mid,x_Base_Width,x_Top_Width,y_Max] = computeStructureBoundary(msh,propALE)

%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Matthew Keller
%
%% Function documentation
%
%  Returns geometrical 
%
%                 Input :
%                   msh : Nodes and elements in the low order mesh:
%                               .nodes : The nodes in the mesh
%                            .elements : The elements in the mesh
%               propALE : Properties regarding the ALE boundary
%                           .nodes : The sequence of the nodal coordinates
%                                    on the ALE boundary
%                       .fcthandle : Function handle to the computation of
%                                    the ALE motion
%                         propUser : Extra user-defined parameters
%
%                Output :
%            x_Base_Min : The minimum x coordinate of the structure
%                 x_Mid : The mid point of the structure for width changes
%          x_Base_Width : The width of the bottom of the structure for
%                         taper calculations
%           x_Top_Width : The width of the top of the structure for
%                         taper calculations
%                 y_Max : The maximum y coordinate of the structure
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
% 3. Return the geometry of structure boundary

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
    
    %% 3. Return the building geometry 
    % Y Geometry  
    y_Max = max(y); % Return max y coordinate of structure boundary
    y_max_index = find(y==y_Max); % Find indices for all nodes along the maximum y boundary
                                  % for top width calculation
    % X Geometry
    x_Base_Min = min(x); % Return min x coordinate of structure boundary
    x_Base_Max = max(x); % Return max x coordinate of structure boundary
    x_Base_Width = x_Base_Max - x_Base_Min; % Return width of structure base
    x_Mid = x_Base_Min + (x_Base_Width*0.5); % Return center of building for dx motion
    x_Top_Min = min(x(y_max_index)); % Return min of the upper boundary nodes
    x_Top_Max = max(x(y_max_index)); % Return max of the upper boundary nodes
    x_Top_Width = x_Top_Max - x_Top_Min; % Return width of structure top
end