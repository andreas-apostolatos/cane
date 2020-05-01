function index = ...
    plot_postprocIGAKirchhoffLoveShell2PatchesLinear ...
    (patch1, dHat1, patch2, dHat2, graph, outMsg)
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
% current configuration of an isogeometric Kirchhoff-Love multipatch shell
% given the Control Point displacement. The second window contains the
% visualization of the selected resultant component over the shell's domain
% in the reference configuration. The visualization handles a two patch
% geometry.
%
%         Input :
% patch1,patch2 : The two patches that form the structure
%   dHat1,dHat2 : The Control Point displacement fields for both patches
%         graph : Information on the graphics
%        outMsg : Whether or not to output message on refinement progress
%                 'outputEnabled' : enables output information
%
%      Output :
%       index : The index of the current graph
%               
%               graphics
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
    fprintf('_________________________________________________________________________\n');
    fprintf('#########################################################################\n');
    fprintf('Plotting postprocessing configuration and resultant visualization for the\n');
    fprintf('linear isogeometric multi-patch Kirchhoff-Love shell has been initiated\n\n');
    fprintf('Configuration to be visualized (1st window): ');
    if strcmp(graph.postprocConfig,'reference')
        fprintf('Reference\n');
    elseif strcmp(graph.postprocConfig,'current')
        fprintf('Current\n');
    elseif strcmp(graph.postprocConfig,'referenceCurrent')
        fprintf('Reference and current\n');
    end
    fprintf('Resultant to be visualized (2nd window): ');
    if strcmp(graph.resultant,'displacement')
        if strcmp(graph.component,'x')
            fprintf('Displacement component u_x\n');
        elseif strcmp(graph.component,'y')
            fprintf('Displacement component u_y\n');
        elseif strcmp(graph.component,'z')
            fprintf('Displacement component u_z\n');
        elseif strcmp(graph.component,'2norm')
            fprintf('Displacement magnitude ||u||_2\n');
        end
    elseif strcmp(graph.resultant,'strain')
        if strcmp(graph.component,'1')
            fprintf('Strain component epsilon_11\n');
        elseif strcmp(graph.component,'2')
            fprintf('Strain component epsilon_22\n');
        elseif strcmp(graph.component,'12')
            fprintf('Strain component epsilon_12\n');
        elseif strcmp(graph.component,'1Principal')
            fprintf('1st pricipal strain field \epsilon_1\n');
        elseif strcmp(graph.component,'2Principal')
            fprintf('2nd pricipal strain field \epsilon_2\n');
        end
    elseif strcmp(graph.resultant,'curvature')
        if strcmp(graph.component,'1')
            fprintf('Curvature component kappa_11\n');
        elseif strcmp(graph.component,'2')
            fprintf('Curvature component kappa_22\n');
        elseif strcmp(graph.component,'12')
            fprintf('Curvature component kappa_12\n');
        elseif strcmp(graph.component,'1Principal')
            fprintf('1st pricipal curvature field \kappa_1\n');
        elseif strcmp(graph.component,'2Principal')
            fprintf('2nd pricipal curvature field \kappa_2\n');
        end
    elseif strcmp(graph.resultant,'force')
        if strcmp(graph.component,'1')
            fprintf('Force component n_11\n');
        elseif strcmp(graph.component,'2')
            fprintf('Force component n_22\n');
        elseif strcmp(graph.component,'12')
            fprintf('Force component n_12\n');
        elseif strcmp(graph.component,'1Principal')
            fprintf('1st pricipal force field n_1\n');
        elseif strcmp(graph.component,'2Principal')
            fprintf('2nd pricipal force field n_2\n');
        end
    elseif strcmp(graph.resultant,'moment')
        if strcmp(graph.component,'1')
            fprintf('Moment component m_11\n');
        elseif strcmp(graph.component,'2')
            fprintf('Moment component m_22\n');
        elseif strcmp(graph.component,'12')
            fprintf('Moment component m_12\n');
        elseif strcmp(graph.component,'1Principal')
            fprintf('1st pricipal moment field m_1\n');
        elseif strcmp(graph.component,'2Principal')
            fprintf('2nd pricipal moment field m_2\n');
        end
    elseif strcmp(graph.resultant,'shearForce')
        if strcmp(graph.component,'1')
            fprintf('Shear force component q_1\n');
        elseif strcmp(graph.component,'2')
            fprintf('Shear force component q_2\n');
        end
    else
        fprintf('No resultant has been chosen to be visualized\n');
    end
    fprintf('_________________________________________________________________________\n\n');
    tic;
end

%% 0. Read input

% Initialize handle to the figure
figure(graph.index)

% For patch 1 :
% _____________

p1 = patch1.p;
q1 = patch1.q;
Xi1 = patch1.Xi;
Eta1 = patch1.Eta;
CP1 = patch1.CP;
isNURBS1 = patch1.isNURBS;
parameters1 = patch1.parameters;
rb1 = patch1.rb;
Fl1 = patch1.Fl;

% For patch 2 :
% _____________

p2 = patch2.p;
q2 = patch2.q;
Xi2 = patch2.Xi;
Eta2 = patch2.Eta;
CP2 = patch2.CP;
isNURBS2 = patch2.isNURBS;
parameters2 = patch2.parameters;
rb2 = patch2.rb;
Fl2 = patch2.Fl;

% Grid point number for the plotting of both the B-Spline surface the knots
% as well as the resultant computation over the domain

% For patch 1 :
% _____________

xiGrid1 = 19;
etaGrid1 = 19;

% For patch 2 :
% _____________

xiGrid2 = 19;
etaGrid2 = 19;

%% 1. Plot the first window: Reference and/or current configuration

% Plot the window
subplot(2,1,1);

% For patch 1 :
% _____________

plot_postprocCurrentConfigurationIGAKirchhoffLoveShell(p1,q1,Xi1,Eta1,CP1,isNURBS1,xiGrid1,etaGrid1,rb1,Fl1,dHat1,graph);
hold on;

% For patch 2 :
% _____________

plot_postprocCurrentConfigurationIGAKirchhoffLoveShell(p2,q2,Xi2,Eta2,CP2,isNURBS2,xiGrid2,etaGrid2,rb2,Fl2,dHat2,graph);

% Assign graphic properties and title
% camlight left; lighting phong;
axis equal;
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
if strcmp(graph.postprocConfig,'reference')
    title('Reference configuration');
elseif strcmp(graph.postprocConfig,'current')
    title('Current configuration/linear');
elseif strcmp(graph.postprocConfig,'referenceCurrent')
    title('Reference and current configuration/linear');
end
hold off;

%% 2. Plot the second window: Resultant visualization

% Plot the window
subplot(2,1,2);

% For patch 1 :
% _____________

plot_postprocResultantsIGAKirchhoffLoveShell(p1,q1,Xi1,Eta1,CP1,isNURBS1,parameters1,xiGrid1,etaGrid1,dHat1,graph);
hold on;

% For patch 2 :
% _____________

plot_postprocResultantsIGAKirchhoffLoveShell(p2,q2,Xi2,Eta2,CP2,isNURBS2,parameters2,xiGrid2,etaGrid2,dHat2,graph);

% Assign graphic properties and title
if strcmp(graph.resultant,'displacement')
    if strcmp(graph.component,'x')
        titleString = 'Displacement component d_x';
    elseif strcmp(graph.component,'y')
        titleString = 'Displacement component d_y';
    elseif strcmp(graph.component,'z')
        titleString = 'Displacement component d_z';
    elseif strcmp(graph.component,'2norm')
        titleString = 'Displacement magnitude ||d||_2';
    end
elseif strcmp(graph.resultant,'strain')
    if strcmp(graph.component,'1')
        titleString = 'Strain component \epsilon_{11}';
    elseif strcmp(graph.component,'2')
        titleString = 'Strain component \epsilon_{22}';
    elseif strcmp(graph.component,'12')
        titleString = 'Strain component \epsilon_{12}';
    elseif strcmp(graph.component,'1Principal')
        titleString = '1st pricipal strain field \epsilon_1';
    elseif strcmp(graph.component,'2Principal')
        titleString = '2nd pricipal strain field \epsilon_2';
    end
elseif strcmp(graph.resultant,'curvature')
    if strcmp(graph.component,'1')
        titleString = 'Curvature component \kappa_{11}';
    elseif strcmp(graph.component,'2')
        titleString = 'Curvature component \kappa_{22}';
    elseif strcmp(graph.component,'12')
        titleString = 'Curvature component \kappa_{12}';
    elseif strcmp(graph.component,'1Principal')
        titleString = '1st pricipal curvature field \kappa_1';
    elseif strcmp(graph.component,'2Principal')
        titleString = '2nd pricipal curvature field \kappa_2';
    end
elseif strcmp(graph.resultant,'force')
    if strcmp(graph.component,'1')
        titleString = 'Force component n_{11}';
    elseif strcmp(graph.component,'2')
        titleString = 'Force component n_{22}';
    elseif strcmp(graph.component,'12')
        titleString = 'Force component n_{12}';
    elseif strcmp(graph.component,'1Principal')
        titleString = '1st pricipal force field n_1';
    elseif strcmp(graph.component,'2Principal')
        titleString = '2nd pricipal force field n_2';
    end
elseif strcmp(graph.resultant,'moment')
    if strcmp(graph.component,'1')
        titleString = 'Moments component m_{11}';
    elseif strcmp(graph.component,'2')
        titleString = 'Moments component m_{22}';
    elseif strcmp(graph.component,'12')
        titleString = 'Moments component m_{12}';
    elseif strcmp(graph.component,'1Principal')
        titleString = '1st pricipal moment field m_1';
    elseif strcmp(graph.component,'2Principal')
        titleString = '2nd pricipal moment field m_2';
    end
elseif strcmp(graph.resultant,'shearForce')
    if strcmp(graph.component,'1')
        titleString = 'Shear force component q_1';
    elseif strcmp(graph.component,'2')
        titleString = 'Shear force component q_2';
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
axis equal;
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);

%% 3. Update the graph index
index = graph.index + 1;

%% 4. Appendix
if strcmp(outMsg, 'outputEnabled')
    computationalTime = toc;
    fprintf('Plotting the current configuration took %.2d seconds \n\n', computationalTime);
    fprintf('__________________Plotting Current Configuration Ended___________________\n');
    fprintf('#########################################################################\n\n\n');
end

end
