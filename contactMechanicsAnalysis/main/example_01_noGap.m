%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Lehrstuhl für Statik, Prof. Dr.-Ing. K.-U. Bletzinger                 %
%   _____________________________________________________                 %
%                                                                         %                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Script documentation
%
% Example contact with no gap 
%
%% Preamble
clc;
clear all;

%% Includes
addpath('../nurbs','../nurbs/nurbs_base_vectors','../nurbs/nurbs_basis_functions',...
        '../nurbs/nurbs_refinement','../nurbs/nurbs_geometry','../mesh',...
        '../mesh2D','../nurbs/nurbs_graphics','../nurbs/nurbs_surface',...
        '../algebra','../nurbs/nurbs_curve','../plot','../load',...
        '../supports',...
        '../boundaryConditions','../planeStressAnalysis','../lagrangeMultipliers');
    
addpath('../../FEMPlateInMembraneActionAnalysis/solutionMatricesAndVectors')
addpath('../../FEMPlateInMembraneActionAnalysis/loads');
addpath('../../basisFunctions');

%% NURBS parameters

% Height
h =2;

% Lenght
l = 4;

% tunel height
tunel_height = 1;

% Inclination of the right edge
inclination = 0;

% Polynomial degrees
p=2;
q=1;

% Knot vectors
U=[0 0 0 1/3 1/3 .5 .5 2/3 2/3 1 1 1];
V=[0 0 1 1];

% Control Point coordinates for the plane structure
CP(:,:,1)= [-(l-1) -(l-1)
            -(l-3) -(l-3)
            -(l-3) -(l-3)
            -(l-3) -(l-3)
            0 0
            (l-3) (l-3)
            (l-3) (l-3)
            (l-3) (l-3)
            (l-1) (l-1)-inclination];
CP(:,:,2)= [0 h
            0 h
            0 h
            tunel_height h
            tunel_height h
            tunel_height h
            0 h
            0 h
            0 h];
CP(:,:,3) = [0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0];
w2 = sin(45);
CP(:,:,4) = [1 1;1 1;1 1; w2 1;1 1;w2 1;1 1; 1 1;1 1];

%% Material constants

% Young's modulus
parameters.E = 1e5;

parameters.rho = 1.0;

% Poisson ratio
parameters.nue = 0.3;

% Thickness of the plate
parameters.t = 1;

%% GUI

% Analysis type -> get from parser
analysis.type = 'plainStrain';

% On the graph
graph.index = 1;
% On the geometry visualization
graph.visualization.geometry = 'current';

%  meshAttributes : Mesh options :
%                    1) boundary : 
%                            uv0 : Meshing parametric line (u,v)=(:,0)
%                            u1v : Meshing parametric line (u,v)=(1,:)
%                            uv1 : Meshing parametric line (u,v)=(:,1)
%                            u0v : Meshing parametric line (u,v)=(0,:)
%                    2) hdata :
%                        hmax : Maximum edge size
%                    3) options : 
%                         maxit : Maximum number of iterations (def. 20)
%                          mlim : The convergence tolerance. The maximum 
%                                 percentage change in edge length per 
%                                 iteration must be less than (def. 0.02)
%                         dhmax : The maximum allowable (relative) gradient 
%                                 in the size function (def.0.3, 30.0%)
%                        output : Displays the mesh and the mesh statistics 
%   `                             upon completion ( def. TRUE)
% 
% Meshing over parametric line (u,v)=(:,0)
meshAttributes.boundary.uv0.number_of_segments = 3;
meshAttributes.boundary.uv0.segment = zeros(meshAttributes.boundary.uv0.number_of_segments,2);
meshAttributes.boundary.uv0.segment(1,:) = [0 1/3];
meshAttributes.boundary.uv0.number_of_nodes(1) = 1;
meshAttributes.boundary.uv0.segment(2,:) = [1/3 2/3];
meshAttributes.boundary.uv0.number_of_nodes(2) = 16;
meshAttributes.boundary.uv0.segment(3,:) = [2/3 1];
meshAttributes.boundary.uv0.number_of_nodes(3) =1;

% Meshing over parametric line (u,v)=(1,:)
meshAttributes.boundary.u1v.number_of_segments = 1;
meshAttributes.boundary.u1v.segment = zeros(meshAttributes.boundary.u1v.number_of_segments,2);
meshAttributes.boundary.u1v.segment(1,:) = [0 1];
meshAttributes.boundary.u1v.number_of_nodes(1) = 1;

% Meshing over parametric line (u,v)=(:,1)
meshAttributes.boundary.uv1.number_of_segments = 1;
meshAttributes.boundary.uv1.segment = zeros(meshAttributes.boundary.uv1.number_of_segments,2);
meshAttributes.boundary.uv1.segment(1,:) = [1 0];
meshAttributes.boundary.uv1.number_of_nodes(1) = 1;

% Meshing over parametric line (u,v)=(0,:)
meshAttributes.boundary.u0v.number_of_segments = 1;
meshAttributes.boundary.u0v.segment = zeros(meshAttributes.boundary.u0v.number_of_segments,2);
meshAttributes.boundary.u0v.segment(1,:) = [1 0];
meshAttributes.boundary.u0v.number_of_nodes(1) = 1;

% Attributes of the mesh
% Maximum edge size
meshAttributes.hdata.hmax  = .15;

% Maximum number of iterations for convergence
meshAttributes.options.maxit = 20;

% Maximum allowable relative gradient in the size function
meshAttributes.options.dhmax = .15;

%% Dirichlet (essential) boundary conditions
boundaryConditions.dirichlet.number_of_segments = 1;
boundaryConditions.dirichlet.number_of_segment = zeros(boundaryConditions.dirichlet.number_of_segments,2);

boundaryConditions.dirichlet.segmentu(1,:) = [0 1/3];
boundaryConditions.dirichlet.segmentv(1,:) = 0;
boundaryConditions.dirichlet.dofs(1) = 3;

boundaryConditions.tolerance_search = 1e-8;
boundaryConditions.increment = 1e-3;
boundaryConditions.tolerance_points = 1e-2;

%% Contact boundary Conditions

% Contact BC can be defined in the same manner as Dirichlet BC.

boundaryConditions.contact.number_of_segments = 1;
boundaryConditions.contact.number_of_segment = zeros(boundaryConditions.contact.number_of_segments,2);

boundaryConditions.contact.segmentu(1,:) = [2/3 1];
boundaryConditions.contact.segmentv(1,:) = 0;

%% rigid wall- line  [(x0,y0) ; (x1,y1)]

wall=[-5, 0; 5 ,0];

segmentsPoints(1,:,:)=wall;


%% Neumann (natural) boundary conditions
boundaryConditions.neumann.number_of_segments = 1;
boundaryConditions.neumann.number_of_segment = zeros(boundaryConditions.neumann.number_of_segments ,2);

boundaryConditions.neumann.segmentu(1,:) = [0 1];
boundaryConditions.neumann.segmentv(1,:) = 1;
boundaryConditions.neumann.load = struct([]);
boundaryConditions.neumann.load(1).type = 'surface_load';
boundaryConditions.neumann.load(1).value = [0 -1e4]';

%% Create a mesh
figure(graph.index)
strMsh = meshCADgeometry(p,U,q,V,CP,meshAttributes);










%% Apply boundary conditions

% get this from GiD
% On the Dirichlet boundary conditions
homDBC = getDirichletBoundaryConditions(p,U,q,V,CP,boundaryConditions,strMsh);

% On the Neumann boundary conditions
F = computeLoadVector(p,U,q,V,CP,boundaryConditions,parameters,strMsh);

% Initialize graphics index
%graph.index = 1;

%% Output data to a VTK format
%pathToOutput = '../../outputVTK/FEMPlateInMembraneActionAnalysis/';

% % Quadrature for the stiffness matrix and the load vector of the problem
% % 'default', 'user'
% intLoad.type = 'default';
% intDomain.type = 'default';
% intLoad.noGP = 1;
% intDomain.noGP = 1;

%% Compute the load vector
%t = 0;
%F = computeLoadVctFEMPlateInMembraneAction(strMsh,analysis,NBC,t,intLoad,'outputEnabled');

%% Visualization of the configuration
%graph.index = plot_referenceConfigurationFEMPlateInMembraneAction...
    %(strMsh,analysis,F,homDBC,graph,'outputEnabled');




% On the Contact boundary (either the whole boundary is a contact area or
% with the function getContactNodes selected regions :

% get this from GiD -> global numbering of nodes
%candidateContactNodes=getContactNodesOnBoundary(p,U,q,V,CP,boundaryConditions,mesh);% Select certain area
candidateContactNodes=strMsh.boundaryNodes;% Define the whole boundary as possible contact area









%% Plot the initial configuration
hold on;
index = plotBoundaryConditionsOnMesh(strMsh,homDBC,F,graph);
plotSegments(strMsh,wall); % plot the wall segment
hold off;

%% Solve the system and get the displacement field
ts=cputime;

maxIterForLagrange=100;

[displacement, lagrange] = solveSignoriniLagrange1(strMsh,homDBC,candidateContactNodes,F,segmentsPoints,parameters,analysis,maxIterForLagrange); % Use Algorithm 1 
%[displacement,lagrange] = solveSignoriniLagrange2(mesh,homDBC,candidateContactNodes,F,segmentsPoints,materialProperties,analysis,maxIterForLagrange); % or use Algorithm 2

fprintf('\t Time :   %4.2f \n',cputime-ts);
%% Postprocessing
graph.index = plotCurrentConfigurationAndResultants(strMsh,homDBC,displacement,graph);
plotSegments(strMsh,segmentsPoints);
plotLagrangeMultipliers( strMsh, displacement,lagrange.active_nodes,lagrange.multipliers,'' );% To show vertical bars for the Lagrange multiplier values insert 'v' as last parameter

%% End of the script