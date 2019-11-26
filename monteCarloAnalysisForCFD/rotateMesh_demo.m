% THIS SCRIPT RUNS FOR A WHILE, AND WRITES ~1MB IMAGES TO THE outputVTK FOLDER
% YOU ARE GOING TO NEED THE export_fig ADDON TO RUN THIS SCRIPT
%
% IF MESH CHECKING IS ENABLED, THIS SCRIPT WON'T WORK
% -> COMMENT OUT THE MESH CHECKING IN rotateBodyInMesh.m
%% SETUP
% Get mesh
problemStruct = init('airfoil_test');
mesh = problemStruct.fldMsh;
% Get body nodes
bodyNodeIndices = getBodyNodes(mesh, problemStruct.homDBC);
%% CONFIG, CALC, EXPORT
% Param set for rotation
angle = -20*pi/180;                                                         % Angle by which the body gets rotated
recursion = 12;
exponent = [ .99 .55 ];                                                     % Angle decrease rate for connected nodes

for rec=1:recursion                                                         % Number of layers of nodes that we want to rotate
    % CALC ----------------------------------------------------------------
    problemStruct.fldMsh.nodes =                                        ...
        rotateBodyInMesh(mesh, bodyNodeIndices, angle, exponent, rec);
    % EXPORT --------------------------------------------------------------
    % Plot
    close all                                                               % Close previous figure
    plotProblem(problemStruct)                                              % Plot new mesh (takes a while)
    set(gcf, 'units','normalized','outerposition',[0 0 1 1])                % Full screen
    axis([-1.5, 3, -1.2 ,1.2])                                              % Set area of interest
    % Assemble file name
    filename = [                                                        ...
        'msh',                                                          ...
        num2str(angle*180/pi, 2),                                       ...
        'deg_exp',                                                      ...
        num2str(exponent, 2),                                           ...
        '_rec',                                                         ...
        num2str(rec, 1),                                                ...
        '.png'                                                          ...
    ];
    filename = [cd, '/outputVTK/', filename];                               % Create unique file name
    % Export figure
    export_fig(filename,'-zbuffer','-r600')                                 % Export high resolution figure (takes a while)
    % ---------------------------------------------------------------------
end
%% END
disp('ready')