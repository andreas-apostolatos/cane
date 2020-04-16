function index = plot_temperatureDistribution ... 
    (mesh, dHat, parameters, propGraph, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%                  Marko Leskovar
%
%% Function documentation
%
% Plots the temperature distribution of the 2D heat transfer problem
%
%             Input :
%            fldMsh : Nodes and elements for the fluid mesh
%        parameters : The parameters of the flow (density, viscosity)
%                 t : Real time at which we want to visualize the results
%         propGraph : Structure containing information on the graphics
%            outMsg : Enables outputting information in the command window 
%                     when selected as 'outputEnabled'
%
%% Function layout :
%
% 0. Read input
%
% 1. Initialize the arrays
%
% 2. Loop over all the evaluation points
%
%    2i. Get the current Cartesian location
%
%    2ii. Compute the resultant on the evaluation point
%
%    2iii. Update the Cartesian loaction along the x-direction
%
% 3. Visualize the selected resultant over the domain
%
% 4. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('____________________________________________________\n');
    fprintf('####################################################\n');
    fprintf('\n');
    fprintf('Temperature Distribution of the heat transfer problem\n');
    fprintf('has been initiated\n');
    fprintf('\n');
    fprintf('____________________________________________________\n');
    fprintf('\n');
    tic;
end

%% 0. Read input






% Initialize figure handle
figure(propGraph.index);
set(figure(propGraph.index),'Name','Temperature Distribution');
hold on

%% 6. Plot the reference configuration

% Define colors
colorDomain = [217 218 219]/255;
colorEdge = 'black';

% Plot the mesh
patch('faces',mesh.elements,'vertices',mesh.nodes,'facecolor',colorDomain,'edgecolor',colorEdge);

%% 7. Plot the temperature distribution

trisurf(mesh.nodes(:,1),mesh.nodes(:,2),dHat);











colorbar('Location','EastOutside');
axis tight;
shading interp;


% Graph properties
axis equal;
xlabel('x', 'FontSize', 14);
ylabel('y', 'FontSize', 14);


hold off;

%% Update the graph index
index = propGraph.index + 1;

%% 4. Appendix
if strcmp(outMsg,'outputEnabled')
    computationalTime = toc;
    fprintf('Plotting the temperature field took %.2d seconds  \n\n', computationalTime);
    fprintf('_________Plotting Numerical Solution Ended__________\n');
    fprintf('####################################################\n');
end

end