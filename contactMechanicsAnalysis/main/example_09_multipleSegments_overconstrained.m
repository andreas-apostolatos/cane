%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Lehrstuhl für Statik, Prof. Dr.-Ing. K.-U. Bletzinger                 %
%   _____________________________________________________                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Script documentation
%
% !!! ATTENTION: This problem is over constrained and can NOT be solved !!!
%
% A wedge-shaped structure between two rigid walls is defined
% and solved with a Lagrange-method. A solver for MULTIPLE contact Signorini 
% problems (no friction) is used 
%
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
h =3;

% x0 upper right corner
x0 = 1;

% upper left corner
x1 = -1;

% lower left corner
x2 = -1;

% lower right corner
x3 = 0.2;

% Polynomial degrees
p=1;
q=1;

% Knot vectors
U=[0 0 1 1];
V=[0 0 1 1];

% Control Point coordinates for the plane structure
CP(:,:,1)= [x2  x1
            x3 x0];
CP(:,:,2)= [0 h
            0 h];
CP(:,:,3) = [0 0; 0 0];
w2 = sin(45);
CP(:,:,4) = [1 1;1 1];

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
graph.visualization.geometry = 'reference_and_current';

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
meshAttributes.boundary.uv0.number_of_segments = 1;
meshAttributes.boundary.uv0.segment = zeros(meshAttributes.boundary.uv0.number_of_segments,2);
meshAttributes.boundary.uv0.segment(1,:) = [0 1];
meshAttributes.boundary.uv0.number_of_nodes(1) = 1;

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
meshAttributes.hdata.hmax  = .1;

% Maximum number of iterations for convergence
meshAttributes.options.maxit = 20;

% Maximum allowable relative gradient in the size function
meshAttributes.options.dhmax = .1;

%% Refinement of the initial geometry
% 
% % Degree elevation 
% tp=0;  tq=0;
% [CP,U,V,p,q] = degree_elevate_surface(p,q,U,V,CP,tp,tq);
% % Knot insertion
% n = 1;
% [CP,U,V] = knot_refinement_consistent_surface(p,U,q,V,CP,n,n);

%% Plot reference geometry
% graph.index = plot_geometry(p,q,U,V,CP,graph);

%% Dirichlet (essential) boundary conditions

boundaryConditions.dirichlet.number_of_segments = 1;
boundaryConditions.dirichlet.number_of_segment = zeros(boundaryConditions.dirichlet.number_of_segments,2);

boundaryConditions.dirichlet.segmentu(1,:) = 0;
boundaryConditions.dirichlet.segmentv(1,:) = [0.92 0.95];
boundaryConditions.dirichlet.dofs(1) = 1;
boundaryConditions.dirichlet.isOnU(1)=0;


%boundaryConditions.tolerance_search = 1e-8;
%boundaryConditions.increment = 1e-3;
%boundaryConditions.tolerance_points = 1e-2;

%% Neumann (natural) boundary conditions
boundaryConditions.neumann.number_of_segments = 1;
boundaryConditions.neumann.number_of_segment = zeros(boundaryConditions.neumann.number_of_segments ,2);

boundaryConditions.neumann.segmentu(1,:) = [0 1];
boundaryConditions.neumann.segmentv(1,:) = 1;
boundaryConditions.neumann.load = struct([]);
boundaryConditions.neumann.load(1).type = 'surface_load';
boundaryConditions.neumann.load(1).value = [0 -1e4]';

%% Contact boundary Conditions

% Here we can define contact boundary conditions in the same manner as 
% the Dirichlet BC. 

% ATTENTION for every rigid wall segment a different set of contact- node 
% candidates should be used !!!!!

% Contact Booundary for wall segment 1
boundaryConditions.contact.number_of_segments = 1;
boundaryConditions.contact.number_of_segment = zeros(boundaryConditions.contact.number_of_segments,2);

boundaryConditions.contact.segmentu(1,:) = 0;
boundaryConditions.contact.segmentv(1,:) = [0 1];

boundaryConditions.tolerance_search = 1e-8;
boundaryConditions.increment = 1e-3;
boundaryConditions.tolerance_points = 1e-2;

% Contact Booundary for wall segment 2
boundaryConditions2.contact.number_of_segments = 1;
boundaryConditions2.contact.number_of_segment = zeros(boundaryConditions2.contact.number_of_segments,2);

boundaryConditions2.contact.segmentu(1,:) = 1;
boundaryConditions2.contact.segmentv(1,:) = [0 1];

boundaryConditions2.tolerance_search = 1e-8;
boundaryConditions2.increment = 1e-3;
boundaryConditions2.tolerance_points = 1e-2;

%% rigid wall- line  [(x0,y0) ; (x1,y1)]


wall1=[x1-(x2-x1)/4,h+h/4;x2+(x2-x1)/1.5,-h/1.5];
wall2=[x3-(x0-x3)/1.5,-h/1.5; x0+(x0-x3)/4,h+h/4];

segmentsPoints(1,:,:)=wall1;
segmentsPoints(2,:,:)=wall2;


%% Create a mesh
figure(graph.index)
mesh = meshCADgeometry(p,U,q,V,CP,meshAttributes);

%% Apply boundary conditions

% On the Dirichlet boundary conditions
rb = getDirichletBoundaryConditions(p,U,q,V,CP,boundaryConditions,mesh);

% On the Neumann boundary conditions
F_global = computeLoadVector(p,U,q,V,CP,boundaryConditions,materialProperties,mesh);

% Define the specific contact boundary for a certain segment
cn1=getContactNodes(p,U,q,V,CP,boundaryConditions,mesh);
cn2=getContactNodes(p,U,q,V,CP,boundaryConditions2,mesh);

% Define the structrure array
cn(1)=struct('indices',cn1);
cn(2)=struct('indices',cn2);

%% Plot the initial configuration
hold on;
index = plotBoundaryConditionsOnMesh(mesh,rb,F_global,graph);
plotSegments(mesh,wall1);
plotSegments(mesh,wall2);
hold off;


%% Solve the system and get the displacement field
ts=cputime;

maxIterForLagrange=5;

[displacement, lagrange] = multSolveSignoriniLagrange1(mesh,rb,cn,F_global,segmentsPoints,materialProperties,analysis,maxIterForLagrange);
%[displacement,lagrange] = multSolveSignoriniLagrange2(mesh,rb,cn,F_global,segmentsPoints,materialProperties,analysis,maxIterForLagrange);

inclination_wall1=h/(x2-x1)

fprintf('\t Time :   %4.2f \n',cputime-ts);

%% Postprocessing
graph.index = plotCurrentConfigurationAndResultants(mesh,rb,displacement,graph);
plotSegments(mesh,segmentsPoints);
plotLagrangeMultipliers( mesh, displacement,lagrange.active_nodes,lagrange.multipliers,'' );

%% End of the scrpit
