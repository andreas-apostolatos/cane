%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Lehrstuhl für Statik, Prof. Dr.-Ing. K.-U. Bletzinger                 %
%   _____________________________________________________                 %
%                                                                         %
%   Authors                                                               %
%   _______                                                               %
%                                                                         %
%   Dr.-Ing. Roland Wüchner                                               %
%   Dipl.-Math. Andreas Apostolatos (andreas.apostolatos@tum.de)          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function index = plotBasisFunctionsTriangularMesh2D(mesh,graph)
%% Function documentation
%
% Plots the basis functions corresponding to the discretization of the 2D
% domain using the linear triangles
%
%   Input :
%    mesh : Structure containing the nodes and the elements of the
%           triangular mesh
%
%  Output :
%   index : The index of the current graph
%
% Function layout :
%
% 1. Graphical output of the mesh
%
% 2. Loop over all elements and plot their basis functions
%
%% Function main body

%% 1. Graphical output of the mesh
figure(graph.index)
axis equal;
cla,patch('vertices',mesh.nodes,'faces',mesh.elements,'edgecol','k','facecol',[.8,.9,1]);

% hold on until to plot all basis functions on the mesh
hold on;

%% 2. Loop over all elements and plot their basis functions

% Grid at each triangular element
grid_xi = 10;
grid_eta = 10;

% step size to both directions
step_xi = 1/(grid_xi+1);
step_eta = 1/(grid_eta+1);

for c_element=1:size(mesh.elements(:,1))

    % Compute the basis functions evaluated at the grid points
    P = zeros(grid_xi+2,grid_eta+2,3);
    N_evaluation = zeros(grid_xi+2,grid_eta+2,3);

    % Initialize local coordinates
    xi = 0;
    eta = 0;
    
    % get element from the mesh
    element = mesh.elements(c_element,:);
    
    % get coordinates of the element vertices
    node_i = element(1,1);
    node_j = element(1,2);
    node_k = element(1,3);
    
    % The vertices of the current triangle
    Pi = mesh.nodes(node_i,:);
    Pj = mesh.nodes(node_j,:);
    Pk = mesh.nodes(node_k,:);

    % The moving vertices
    Pi_eta = Pi;
    Pj_eta = Pj;
    
    % Loop over all the sampling points
    for j=1:grid_eta+2
        for i=1:grid_xi+2
            % Get the point in the interior of the line defined from Pi and Pj
            P(i,j,:) = xi*Pi_eta+(1-xi)*Pj_eta;
        
            % Evaluate the basis functions at this point. Elevate their values
            % by the physical elevation of the mesh in the z-direction
            N_evaluation(i,j,:) = P(i,j,3) + computeBasisFunctionsTriangle2D(Pi,Pj,Pk,P(i,j,1),P(i,j,2));
        
            % Update local coordinate by the step size
            xi = xi + step_xi;
        end
        % Reset xi local coordinate
        xi =0;
    
    	% Update eta local coordinate
        eta = eta + step_eta;
    
        % Update the moving vertices over the triangle edges
        Pi_eta = eta*Pk+(1-eta)*Pi;
        Pj_eta = eta*Pk+(1-eta)*Pj;
    end


    % Graphical output  of the basis functions
    figure(graph.index)
    surf(P(:,:,1),P(:,:,2),N_evaluation(:,:,1),'FaceColor','red','EdgeColor','none');
    hold on;
    surf(P(:,:,1),P(:,:,2),N_evaluation(:,:,2),'FaceColor','blue','EdgeColor','none');
    hold on;
    surf(P(:,:,1),P(:,:,2),N_evaluation(:,:,3),'FaceColor','black','EdgeColor','none');
    hold on;

end

% graph properties
drawnow
view (2);
grid on;
hold off;

% Update graph index
index = graph.index + 1;

end

