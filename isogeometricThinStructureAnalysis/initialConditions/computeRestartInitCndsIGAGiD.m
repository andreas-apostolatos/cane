function [dHatHistoryRestart,dHatDotRestart,dHatDDotRestart,noTimeStepRestart,propTransientAnalysis] = ...
    computeRestartInitCndsIGAGiD...
    (BSplinePatches,noDOFs,propTransientAnalysis,caseName,pathToOutput,tab,outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the displacement, the velocity and the acceleration fields
% corresponding to the restart time step retrieved from an output file.
%
%                 Input :
%        BSplinePatches : Cell array containing the B-Spline patches which
%                         are forming the multipatch geometry
%                noDOFs : Number of DOFs for the mulitpatch system 
%                         including also Lagrange Multipliers DOFs if they
%                         are employed
% propTransientAnalysis : Structure containing information on the transient
%                         analysis
%                               .method : 'explicitEuler', 'bossak' etc.
%                               .TStart : Start time of the simulation
%                                 .TEnd : End time of the simulation
%                          .noTimeSteps : Number of time steps
%                    .noTimeStepRestart : Number of the time step from
%                                         which to restart
%              caseName : Dummy variable for this function
%          pathToOutput : Dummy variable for this function
%                   tab : Tabulation for outputting message in the command
%                         window
%                outMsg : Enables outputting information on the command
%                         window when chosen as 'outputEnabled'
%
%                Output :
%    dHatHistoryRestart : The history of the displacement field until the
%                         restart time step
%        dHatDotRestart : The velocity at the restart time step
%       dHatDDotRestart : The acceleration at the restart time step
%     noTimeStepRestart : The number of the starting time step
% propTransientAnalysis : The updated transient analysis properties
%                         structure
%
% Function layout :
%
% 0. Check input
%
% 1. Compute the initial conditions
%
% 2. Loop over all time steps
% ->
%    2i. Make the counter of the time step string
%
%   2ii. Save the discrete displacement, velocity and acceleration field
%
%  2iii. Get the solution vector
%
%   2iv. Save the displacement of the current time step
%
%    2v. Update the time derivatives of the field
% <-
%
%% Function main body

%% 0. Check input

% Time step information on the output file
string1 = 'Result "displacement" "time_step"';
string2 = 'Vector OnNurbsSurface';
string3 = 'Values'; 
string4 = 'End Values';
string5 = 'Vector OnNurbsSurface';
string6 = '"time_step"';

% Number of patches
noPatches = length(BSplinePatches);

% Check if output file is empty
outputFile = fileread([pathToOutput caseName '.gid' '/' caseName '.post.res']); 
if isempty(outputFile)
    error('The output file to retrieve restart data is empty');
end

% Find at which time step to restart
textCut = regexp(outputFile,string4,'split');
textCut2 = regexp(textCut{end - 1},string5,'split');
textCut3 = regexp(textCut2{1},string6,'split');
noTimeStepStr = textCut3{2};
noTimeStepRestart = str2double(noTimeStepStr);
if strcmp(outMsg,'outputEnabled')
    fprintf(strcat(tab,'Restarting at time step %d\n'),noTimeStepRestart);
end

% Initialize the output array
dHatHistoryRestart = zeros(noDOFs,noTimeStepRestart + 1);

%% 1. Compute the initial conditions
warning('The restart option takes into account null initial displacement, velocity and acceleration');
[dHat,dHatDotRestart,dHatDDotRestart,~] = computeNullInitCndsIGAThinStructure...
    (BSplinePatches,noDOFs,propTransientAnalysis,caseName,pathToOutput,tab,outMsg);

%% 2. Loop over all time steps
for iTimeStep = 0:noTimeStepRestart
    % iTimeStep
    %% 2i. Make the counter of the time step string
    iTimeStepStr = num2str(iTimeStep);
    
    %% 2ii. Save the discrete displacement, velocity and acceleration field
    dHatSaved = dHat;
    dHatDotSaved = dHatDotRestart;
    dHatDDotSaved = dHatDDotRestart;
    
    %% 2iii. Get the solution vector
    
    % Split the text accrording to the time step number
    timeStepRestart = regexp(outputFile,[string1 ' ' iTimeStepStr ' ' string2],'split');
    timeStepRestart = regexp(timeStepRestart{2},string3,'split');
    timeStepRestart = timeStepRestart{2};
    timeStepRestart = regexp(timeStepRestart,string4,'split');
    timeStepRestart = timeStepRestart{1};

    % Get the displacement vector containing the patch numbering
    dHat = textscan(timeStepRestart,'%f','delimiter',' ','MultipleDelimsAsOne', 1);
    dHat = dHat{1};

    % Erase the patch numbering from the displacement vector
    dHat(1,:) = [];
    noDOFs = 0;
    for iPatches = 1:noPatches - 1
        noDOFsPatch = BSplinePatches{iPatches}.noDOFs;
        noDOFs = noDOFs + noDOFsPatch;
        dHat(noDOFs + 1,:) = [];
    end
    
    %% 2iv. Save the displacement of the current time step
    dHatHistoryRestart(:,iTimeStep + 1) = dHat;
    
    %% 2v. Update the time derivatives of the field
    if isa(propTransientAnalysis.computeUpdatedVct,'function_handle')
        [dHatDotRestart,dHatDDotRestart] = propTransientAnalysis.computeUpdatedVct ...
            (dHat,dHatSaved,dHatDotSaved,dHatDDotSaved,propTransientAnalysis);
    else
        error('Function handle propTransientAnalysis.computeUpdatedVct undefined');
    end
end

end
