function index = plot_postprocIGAMembraneMultipatchesNLinear ...
    (BSplinePatches, dHat, graph, outMsg)
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
% current configuration of an isogeometric multipatch membrane given the 
% Control Point displacement. The second window contains the visualization 
% of the selected resultant component over the shell's domain in the 
% reference configuration. The visualization handles a two patch geometry. 
% The visualization assumes a nonlinear strain measure and accounts for an 
% arbitray number of membrane multipatches of the coupled system.
%
%          Input :
% BSplinePatches : Structure array containing all the patches of the coupled
%                  system
%           dHat : Structure array containing the displacement field for
%                  each patch of the coupled system
%          graph : Information on the graphics
%         outMsg : Whether or not to output message on refinement progress
%                  'outputEnabled' : enables output information
%
%       Output :
%        index : The index of the current graph
%               
%                graphics
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
    fprintf('__________________________________________________________________________\n');
    fprintf('##########################################################################\n');
    fprintf('Plotting postprocessing configuration and resultant visualization for the\n');
    fprintf('nonlinear isogeometric multi-patch membrane has been initiated\n\n');
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
    elseif strcmp(graph.resultant,'force')
        if strcmp(graph.component,'1')
            fprintf('Force component n^11\n');
        elseif strcmp(graph.component,'2')
            fprintf('Force component n^22\n');
        elseif strcmp(graph.component,'12')
            fprintf('Force component n^12\n');
        elseif strcmp(graph.component,'1Principal')
            fprintf('1st pricipal force field n^1\n');
        elseif strcmp(graph.component,'2Principal')
            fprintf('2nd pricipal force field n^2\n');
        end
    else
        error('No resultant has been chosen to be visualized\n');
    end
    fprintf('__________________________________________________________________________\n\n');
    tic;
end

%% 0. Read input

% Initialize handle to the figure
figure(graph.index)

% Number of patches in the coupled system
noPatches = length(BSplinePatches);

% Check input
for iPatches = 1:noPatches
    if isfield(BSplinePatches{iPatches}.weakDBC,'computeTangMtxResVct') && ...
            (strcmp(BSplinePatches{iPatches}.weakDBC.method,'penalty') || ...
            strcmp(BSplinePatches{iPatches}.weakDBC.method,'lagrangeMultipliers'))
        error('No field weakDBC.%s needs to be defined for the method %s corresponding to the application of weak Dirichlet boundary conditions', ...
            BSplinePatches{iPatches}.weakDBC.computeTangMtxResVct,BSplinePatches{iPatches}.weakDBC.method);
    end
end

% Grid point number for the plotting of both the B-Spline surface the knots
% as well as the resultant computation over the domain
xiGrid = 39;
etaGrid = 39;

% Create an element freedom tables for each patch
for i = 1:noPatches
    if i == 1
        BSplinePatches{i}.EFTPatches = 1:3*BSplinePatches{i}.noCPs;
    else
        BSplinePatches{i}.EFTPatches = ...
            BSplinePatches{i-1}.EFTPatches(length(BSplinePatches{i-1}.EFTPatches)) + ...
            1:BSplinePatches{i-1}.EFTPatches(length(BSplinePatches{i-1}.EFTPatches)) + ...
            3*BSplinePatches{i}.noCPs;
    end
end

%% 1. Plot the first window: Reference and/or current configuration

% Plot the window
subplot(2,1,1);

% Loop over all the patches in the coupled system
for iPatches = 1:noPatches
    if ~exist('BSplinePatches{counterPatches}.FGamma','var')
        BSplinePatches{iPatches}.FGamma = ...
            zeros(3*BSplinePatches{iPatches}.noCPs,1);
    end
    plot_postprocCurrentConfigurationIGAThinStructure(BSplinePatches{iPatches}.p,BSplinePatches{iPatches}.q,BSplinePatches{iPatches}.Xi,...
        BSplinePatches{iPatches}.Eta,BSplinePatches{iPatches}.CP,BSplinePatches{iPatches}.isNURBS,xiGrid,etaGrid,BSplinePatches{iPatches}.homDOFs,...
        BSplinePatches{iPatches}.FGamma,dHat(BSplinePatches{iPatches}.EFTPatches),graph);
    hold on;
end

% Assign graphic properties and title
% camlight left;
lighting phong;
axis equal;
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
zlabel('z','FontSize',14);
if strcmp(graph.postprocConfig,'reference')
    title('Reference configuration');
elseif strcmp(graph.postprocConfig,'current')
    title('Current configuration/linear');
elseif strcmp(graph.postprocConfig,'referenceCurrent')
    title('Reference and current configuration/linear');
end
hold off;

% axis off;
title '';

%% 2. Plot the second window: Resultant visualization

% Plot the window
subplot(2,1,2);

% Loop over all the patches in the coupled system
for iPatches = 1:noPatches
    plot_postprocResultantsIGAMembraneNLinear(BSplinePatches{iPatches}.p,BSplinePatches{iPatches}.q,BSplinePatches{iPatches}.Xi,...
        BSplinePatches{iPatches}.Eta,BSplinePatches{iPatches}.CP,BSplinePatches{iPatches}.isNURBS,BSplinePatches{iPatches}.parameters,...
        BSplinePatches{iPatches}.DOFNumbering,xiGrid,etaGrid,dHat(BSplinePatches{iPatches}.EFTPatches),graph);
    hold on;
end

% Assign graphic properties and title
if strcmp(graph.resultant,'displacement')
    if strcmp(graph.component,'x')
        titleString = 'Displacement component d_x';
    elseif strcmp(graph.component,'y')
        titleString = 'Displacement component d_y';
    elseif strcmp(graph.component,'z')
        titleString = 'Displacement component d_z';
    elseif strcmp(graph.component,'2norm')
        titleString = 'Displacement magnitude || d ||_2';
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
elseif strcmp(graph.resultant,'force')
    if strcmp(graph.component,'1')
        titleString = 'Force component n^{11}';
    elseif strcmp(graph.component,'2')
        titleString = 'Force component n^{22}';
    elseif strcmp(graph.component,'12')
        titleString = 'Force component n^{12}';
    elseif strcmp(graph.component,'1Principal')
        titleString = '1st pricipal force field n^1';
    elseif strcmp(graph.component,'2Principal')
        titleString = '2nd pricipal force field n^2';
    end
end
title(titleString);
shading interp;
colormap('jet');

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
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    fprintf('Plotting the current configuration took %.2d seconds \n\n',computationalTime);
    fprintf('__________________Plotting Current Configuration Ended____________________\n');
    fprintf('##########################################################################\n\n\n');
end

end
