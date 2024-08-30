function index = plot_postprocIGAPlateInMembraneAction ...
    (p, q, Xi, Eta, CP, isNURBS, homDOFs, propParameters, Fl, dHat, ...
    propGraph, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Plots two windows in one: The first contains the reference and/or the
% current configuration of an isogeometric plate in membrane action given
% the Control Point displacement. The second window contains the
% visualization of the selected resultant component over the plate's domain
% in the reference configuration.
%
%          Input :
%            p,q : Polynomial degrees
%         Xi,Eta : Knot vectors in xi,eta-direction
%             CP : Control point coordinates and weights of the undeformed 
%                  plate
%        isNURBS : Flag on whether the geometrical basis is NURBS or 
%                  B-Spline
%        homDOFs : Global numbering of the DOFs where homogeneous Dirichlet
%                  boundary conditions are applied
% propParameters : Structure containing information on the material
%                  parameters,
%                   .E : Young's modulus
%                 .nue : Poisson's ratio
%                    t : thickness
%             Fl : The applied load vector
%           dHat : The displacement field of the control points
%      propGraph : Structure containing information on the figures,
%                      .index : Index of the current figure
%         outMsg : Whether or not to output message on refinement progress
%                  'outputEnabled' : enables output information
%
%         Output :
%          index : The index of the current graph
%               
%                  graphics
%
% Function layout :
%
% 0. Read input
%
% 1. Plot the first window: Reference and/or current configuration
%
% 2. Plot the second window: Resultant visualization
%
% 3. Update the graph index
%
% 4. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('__________________________________________________________________\n');
    fprintf('##################################################################\n');
    fprintf('Plotting postprocessing configuration and resultant visualization\n');
    fprintf('for the isogeometric plate in mebrane action has been initiated\n\n');
    fprintf('Configuration to be visualized (1st window): ');
    if strcmp(propGraph.postprocConfig, 'reference')
        fprintf('Reference\n');
    elseif strcmp(propGraph.postprocConfig, 'current')
        fprintf('Current\n');
    elseif strcmp(propGraph.postprocConfig, 'referenceCurrent')
        fprintf('Reference and current\n');
    end
    fprintf('Resultant to be visualized (2nd window): ');
    if strcmp(propGraph.resultant,'displacement')
        if strcmp(propGraph.component, 'x')
            fprintf('Displacement component u_x\n');
        elseif strcmp(propGraph.component, 'y')
            fprintf('Displacement component u_y\n');
        elseif strcmp(propGraph.component, '2norm')
            fprintf('Displacement magnitude ||u||_2\n');
        end
    elseif strcmp(propGraph.resultant, 'strain')
        if strcmp(propGraph.component, 'x')
            fprintf('Strain component epsilon_xx\n');
        elseif strcmp(propGraph.component, 'y')
            fprintf('Strain component epsilon_yy\n');
        elseif strcmp(propGraph.component, 'xy')
            fprintf('Strain component epsilon_xy\n');
        end
    elseif strcmp(propGraph.resultant, 'stress')
        if strcmp(propGraph.component, 'x')
            fprintf('Stress component sigma_xx\n');
        elseif strcmp(propGraph.component, 'y')
            fprintf('Stress component sigma_yy\n');
        elseif strcmp(propGraph.component, 'xy')
            fprintf('Stress component sigma_xy\n');
        end
    end
    fprintf('__________________________________________________________________\n\n');
    tic;
end

%% 0. Read input

% Initialize handle to the figure
figure(propGraph.index)

% Grid point number for the plotting of both the B-Spline surface the knots
% as well as the resultant computation over the domain
xiGrid = 49;
etaGrid = 49;

%% 1. Plot the first window: Reference and/or current configuration

% Plot the window
subplot(2, 1, 1);
plot_postprocCurrentConfigurationIGAPlateInMembraneAction ...
    (p, q, Xi, Eta, CP, isNURBS, xiGrid, etaGrid, homDOFs, Fl, ...
    dHat, propGraph);

% Assign graphic properties and title
camlight left; 
lighting phong;
view(2);
axis equal;
xlabel('x', 'FontSize', 14);
ylabel('y', 'FontSize', 14);
if strcmp(propGraph.postprocConfig, 'reference')
    title('Reference configuration');
elseif strcmp(propGraph.postprocConfig, 'current')
    title('Current configuration');
elseif strcmp(propGraph.postprocConfig, 'referenceCurrent')
    title('Reference and current configuration');
end
hold off;

%% 2. Plot the second window: Resultant visualization

% Plot the window
subplot(2, 1, 2);
plot_postprocResultantsIGAPlateInMembraneAction ...
    (p, q, Xi, Eta, CP, isNURBS, propParameters, xiGrid, etaGrid, ...
    dHat, propGraph);

% Assign graphic properties and title
if strcmp(propGraph.resultant, 'displacement')
    if strcmp(propGraph.component, 'x')
        titleString = 'Displacement component d_x';
    elseif strcmp(propGraph.component, 'y')
        titleString = 'Displacement component d_y';
    elseif strcmp(propGraph.component, '2norm')
        titleString = 'Displacement magnitude ||d||_2';
    end
elseif strcmp(propGraph.resultant, 'strain')
    if strcmp(propGraph.component, 'x')
        titleString = 'Strain component \epsilon_{xx}';
    elseif strcmp(propGraph.component, 'y')
        titleString = 'Strain component \epsilon_{yy}';
    elseif strcmp(propGraph.component, 'xy')
        titleString = 'Strain component \epsilon_{xy}';
    elseif strcmp(propGraph.component, '1Principal')
        titleString = '1st pricipal strain field \epsilon_1';
    elseif strcmp(propGraph.component, '2Principal')
        titleString = '2nd pricipal strain field \epsilon_2';
    end
elseif strcmp(propGraph.resultant, 'stress')
    if strcmp(propGraph.component, 'x')
        titleString = 'Stress component \sigma_{xx}';
    elseif strcmp(propGraph.component, 'y')
        titleString = 'Stress component \sigma_{yy}';
    elseif strcmp(propGraph.component, 'xy')
        titleString = 'Stress component \sigma_{xy}';
    elseif strcmp(propGraph.component, '1Principal')
        titleString = '1st pricipal stress field \sigma_1';
    elseif strcmp(propGraph.component, '2Principal')
        titleString = '2nd pricipal stress field \sigma_2';
    end
end
title(titleString);

shading interp;
colormap('default');

% invert default colormap => red = negativ, blue = positive
% COL = colormap;
% invCOL(:,1) = COL(:,3);
% invCOL(:,2) = COL(:,2);
% invCOL(:,3) = COL(:,1);
% colormap(invCOL);

% make colormap symmetric
% colim = caxis;
% caxis([-max(abs(colim)) max(abs(colim))]);
colorbar;

view(2);
axis equal;
xlabel('x', 'FontSize', 14);
ylabel('y', 'FontSize', 14);

%% 3. Update the graph index
index = propGraph.index + 1;

%% 4. Appendix
if strcmp(outMsg, 'outputEnabled')
    computationalTime = toc;
    fprintf('Plotting the current configuration took %.2d seconds \n\n', computationalTime);
    fprintf('______________Plotting Current Configuration Ended________________\n');
    fprintf('##################################################################\n\n\n');
end

end