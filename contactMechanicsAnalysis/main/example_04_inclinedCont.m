%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Lehrstuhl für Statik, Prof. Dr.-Ing. K.-U. Bletzinger                 %
%   _____________________________________________________                 %
%                                                                         %                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Script documentation
%
% Example: Inclined segment
%
%% Preamble
clc;
clear all;

%% Includes
addpath('../nurbs','../nurbs/nurbs_base_vectors','../nurbs/nurbs_basis_functions',...
        '../nurbs/nurbs_refinement','../nurbs/nurbs_geometry','../mesh',...
        '../mesh2D','../nurbs/nurbs_graphics','../nurbs/nurbs_surface',...
        '../algebra','../nurbs/nurbs_curve','../basisFunctions','../plot','../load',...
        '../supports','../stiffnessMatrices','../plane_stress_analysis',...
        '../boundaryConditions','../planeStressAnalysis','../lagrangeMultipliers');

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
materialProperties.E = 1e5;

% Poisson ratio
materialProperties.nu = 0.3;

% Thickness of the plate
materialProperties.t = 1;

%% GUI

% Analysis type
analysis.dimension = '2d';
analysis.dofs = 'displacements';
analysis.physics = 'plain_strain';
analysis.type = 'linear';

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

wall=[0, -1; 4 ,-0.5];

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
mesh = meshCADgeometry(p,U,q,V,CP,meshAttributes);

%% Apply boundary conditions

% On the Dirichlet boundary conditions
rb = getDirichletBoundaryConditions(p,U,q,V,CP,boundaryConditions,mesh);

% On the Neumann boundary conditions
F_global = computeLoadVector(p,U,q,V,CP,boundaryConditions,materialProperties,mesh);

% On the Contact boundary (either the whole boundary is a contact area or
% with the function getContactNodes selected regions :

%cn=getContactNodesOnBoundary(p,U,q,V,CP,boundaryConditions,mesh);% Select certain area
cn=mesh.boundaryNodes;% Define the whole boundary as possible contact area

%% Plot the initial configuration
hold on;
index = plotBoundaryConditionsOnMesh(mesh,rb,F_global,graph);
plotSegments(mesh,wall); %plot the wall segment
hold off;

%% Solve the system and get the displacement field
ts=cputime;

maxIterForLagrange=100;

[displacement, lagrange] = solveSignoriniLagrange1(mesh,rb,cn,F_global,segmentsPoints,materialProperties,analysis,maxIterForLagrange); % Use Algorithm 1 
%[displacement,lagrange] = solveSignoriniLagrange2(mesh,rb,cn,F_global,segmentsPoints,materialProperties,analysis,maxIterForLagrange); % or use Algorithm 2

fprintf('\t Time :   %4.2f \n',cputime-ts);
%% Postprocessing
graph.index = plotCurrentConfigurationAndResultants(mesh,rb,displacement,graph);
plotSegments(mesh,segmentsPoints);
plotLagrangeMultipliers( mesh, displacement,lagrange.active_nodes,lagrange.multipliers,'' );% To show vertical bars for the Lagrange multiplier values insert 'v' as last parameter

%% End of the scrpit