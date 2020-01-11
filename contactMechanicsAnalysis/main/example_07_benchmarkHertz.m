%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Script documentation
%
%Contact Problem: a half cylinder with a line Force. The problem is solved
% analytically by Hertz-Theory
%
%% Preamble

clear all;
clc;
ts=cputime;
%% Includes
addpath('../nurbs','../nurbs/nurbs_base_vectors','../nurbs/nurbs_basis_functions',...
        '../nurbs/nurbs_refinement','../nurbs/nurbs_geometry','../mesh',...
        '../mesh2D','../nurbs/nurbs_graphics','../nurbs/nurbs_surface',...
        '../algebra','../nurbs/nurbs_curve','../basisFunctions','../plot','../load',...
        '../supports','../stiffnessMatrices','../plane_stress_analysis',...
        '../boundaryConditions','../planeStressAnalysis','../lagrangeMultipliers');

%% NURBS parameters and point force

% Radius
r =5;

% line Force on Top
F = 1e4;


% Polynomial degrees
p=2;
q=1;

% Knot vectors
U=[0 0 0 0.5 0.5 1 1 1];
V=[0 0 1 1];

 
% Control Point coordinates for the plane structure
CP(:,:,1)= [-r -r   
            -r -r
            0 0
            r r
            r r];
CP(:,:,2)= [0 0
            0 r
            0 r
            0 r
            0 0];
CP(:,:,3) = [0 0; 0 0;0 0;0 0;0 0];
w2 = sin(pi/4);
CP(:,:,4) = [1 1;w2 w2;1 1;w2 w2;1 1];






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
meshAttributes.boundary.uv1.number_of_nodes(1) = 30;

% Meshing over parametric line (u,v)=(0,:)
meshAttributes.boundary.u0v.number_of_segments = 1;
meshAttributes.boundary.u0v.segment = zeros(meshAttributes.boundary.u0v.number_of_segments,2);
meshAttributes.boundary.u0v.segment(1,:) = [1 0];
meshAttributes.boundary.u0v.number_of_nodes(1) = 15;

% Attributes of the mes
% Maximum edge size
meshAttributes.hdata.hmax  = 0.3;

% Maximum number of iterations for convergence
meshAttributes.options.maxit = 20;

% Maximum allowable relative gradient in the size function
meshAttributes.options.dhmax = .3;

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

boundaryConditions.dirichlet.isOnU(1) = 1;
boundaryConditions.dirichlet.segmentu(1,:) = [0 1];
boundaryConditions.dirichlet.segmentv(1,:) = 0;
boundaryConditions.dirichlet.dofs(1) = 2;

boundaryConditions.tolerance_search = 1e-8;
boundaryConditions.increment = 1e-4;
boundaryConditions.tolerance_points = 1e-1;


%% Contact boundary Conditions
boundaryConditions.contact.number_of_segments = 1;
boundaryConditions.contact.number_of_segment = zeros(boundaryConditions.contact.number_of_segments,2);

boundaryConditions.contact.segmentu(1,:) = [0.75 1];
boundaryConditions.contact.segmentv(1,:) = 1;

%% rigid wall- line  [(x0,y0) ; (x1,y1)]

wall=[r, -r/2; r ,r];

segmentsPoints(1,:,:)=wall;

%% Neumann (natural) boundary conditions
boundaryConditions.neumann.number_of_segments = 1;
boundaryConditions.neumann.number_of_segment = zeros(boundaryConditions.neumann.number_of_segments ,2);


boundaryConditions.neumann.segmentu(1) = 0;
boundaryConditions.neumann.segmentv(1) = 0;

boundaryConditions.neumann.isOnU(1) = 1;
boundaryConditions.neumann.load = struct([]);
boundaryConditions.neumann.load(1).type = 'line_load';   % ATTENTION a line load was defined
boundaryConditions.neumann.load(1).value = [F 0]';

%% Create a mesh
figure(graph.index)
mesh = meshCADgeometry(p,U,q,V,CP,meshAttributes);

%% Apply boundary conditions

% On the Dirichlet boundary conditions
rb = getDirichletBoundaryConditions(p,U,q,V,CP,boundaryConditions,mesh);

% On the Neumann boundary conditions
%gl=getDirichletBoundaryConditions(p,U,q,V,CP,loadBoundaryConditions,mesh);

F_global = computeLoadVector(p,U,q,V,CP,boundaryConditions,materialProperties,mesh);

% On the Contact boundary
%cn=getContactNodesOnBoundary(p,U,q,V,CP,boundaryConditions,mesh);
cn=mesh.boundaryNodes;
%% Plot the initial configuration
hold on;
index = plotBoundaryConditionsOnMesh(mesh,rb,F_global,graph);
plotSegments(mesh,wall);
hold off;



%% Solve the system and get the displacement field
ts=cputime;

maxIterForLagrange=100;

[displacement, lagrange] = solveSignoriniLagrange1(mesh,rb,cn,F_global,segmentsPoints,materialProperties,analysis,maxIterForLagrange);
%[displacement,lagrange] = solveSignoriniLagrange2(mesh,rb,cn,F_global,segmentsPoints,materialProperties,analysis,maxIterForLagrange);

fprintf('\t Time :   %4.2f \n',cputime-ts);
%% Postprocessing
graph.index = plotCurrentConfigurationAndResultants(mesh,rb,displacement,graph);
plotSegments(mesh,segmentsPoints);
plotLagrangeMultipliers( mesh, displacement,lagrange.active_nodes,lagrange.multipliers,'h' );% To show vertical bars for the Lagrange multiplier values insert 'v' as last parameter

%% Get the length of the contact area and the reaction force on the contact
[contact_length,contact_force,max_contact_press] = contactLengthForcePressure(mesh,displacement,cn,lagrange.multipliers,materialProperties);

hertz_contact_length=sqrt(4*(2*F)*r*((1-materialProperties.nu^2)/materialProperties.E)/(pi*materialProperties.t));
hertz_pressure=2*(2*F)/(materialProperties.t*pi*hertz_contact_length);


fprintf('\n \n \n \t                    The COMPUTED contact length is :   %f \n',contact_length);
fprintf('\t   The contact length according to HERTZ theory is :   %f \n\n',hertz_contact_length);

fprintf('\t                   The maximal COMPUTED pressure is:   %f \n',max_contact_press);
fprintf('\t The maximal pressure according to HERTZ theory is :   %f \n',hertz_pressure);

fprintf('\t Time :   %4.2f \n',cputime-ts);




%% End of the scrpit

