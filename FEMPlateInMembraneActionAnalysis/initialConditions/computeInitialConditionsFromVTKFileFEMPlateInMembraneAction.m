function [d, dDot, dDDot, numTimeStep] = ...
    computeInitialConditionsFromVTKFileFEMPlateInMembraneAction ...
    (propAnalysis, strMsh, DOF4Output, propParameters, strDynamics, ...
    VTKResultFile, caseName, pathToFile)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the initial displacement field as well as its first and second 
% order time derivative for a CSD simulation which has been stopped, by 
% reading the results from the generated output files defined in 
% VTKResultFile.
%
%             Input :
%      propAnalysis : Structure containing general information on the 
%                     analysis,
%                       .type : Analysis type
%            fldMsh : Nodes and elements of the structural mesh
%        DOF4Output : Array containing the arrangment of the DOFs for
%                     printing them out
%    propParameters : Material parameters (null variable for this function)
%       strDynamics : Transient analysis parameters for the CSD simulation
%     VTKResultFile : The name of the result file in the output folder
%                     where to get the initial conditions for the transient 
%                     simulation
%          caseName : The name of the case
%        pathToFile : Path to the file where to get the initial conditions 
%                     from
%
%            Output :
%                 d : Vector of DOFs containing the initial conditions for
%                     the displacement field
%              dDot : Vector of DOFs containing the initial conditions for 
%                     the velocity field
%             dDDot : Vector of DOFs containing the initial conditions for 
%                     the acceleration field
%       numTimeStep : The time step where the simulation has been paused
%
% Function layout :
%
% 0. Read input
%
% 1. Extract the time step number from the current and the previous time step
%
% 2. Compute the starting time of the simulation
%
% 3. Get the displacement field of the current time step
%
% 4. Get the rates of the displacement field at the current time step
%
%% Function main body

%% 0. Read input

% Number of nodes
numNodes = length(strMsh.nodes(:, 1));

% Number of degees of freedom
if strcmp(propAnalysis.type, 'planeStress') || ...
        strcmp(propAnalysis.type, 'planeStrain')
    numDOFs = 2*numNodes;
else
    error('Analysis type %s is not supported in this function', ...
        propAnalysis.type);
end

% Get the time step when the simulation was stopped
[~, vtkFileName, ~] = fileparts(strcat(pathToFile, VTKResultFile, '.vtk'));

% Initialize output arrays
d = zeros(numDOFs, 1);

%% 1. Extract the time step number from the current and the previous time step
counter = 1;
for i = length(vtkFileName):-1:1
    if ~strcmp(vtkFileName(i), '_')
        timeStepNoString(counter) = vtkFileName(i);
        counter = counter + 1;
    else
        break;
    end
end
timeStepNoString = fliplr(timeStepNoString);
timeStepNo = str2double(timeStepNoString);

%% 2. Compute the starting time of the simulation
numTimeStep = timeStepNo - 1;

%% 3. Get the displacement field of the current time step

% Read the vtk output file for the current time step
vtkFileNameCurrentTimeStep = ...
    fileread(strcat(pathToFile, caseName, '/', caseName, VTKResultFile, '.vtk'));

% Read in the displacement field
block = regexp(vtkFileNameCurrentTimeStep, 'VECTORS displacements double', 'split');
block(1) = [];
out = cell(size(block));
for k = 1:numel(block)
    out{k} = textscan(block{k},'%f %f %f','delimiter',' ','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
velocities = cell2mat(out);

% Assemble the vector of the initial conditions for the unknown
% displacement field
for i = 1:2
    d(DOF4Output(i, :)) = velocities(:, i);
end

%% 4. Get the rates of the displacement field at the current time step

% Get the file name of the VTK result file for the rates of the current
% time step
VTKResultFileRatesCurrentTimeStep = strcat('_dRates_', timeStepNoString);

% Read the file
vtkFileNameRatesCurrentTime = ...
    fileread(strcat(pathToFile, caseName, '/dRates/', caseName, ...
    VTKResultFileRatesCurrentTimeStep, '.txt'));

% Read in the vector of the velocity DOFs
block = regexp(vtkFileNameRatesCurrentTime, 'dDot ', 'split');
block(1) = [];
out = cell(size(block));
for k = 1:numel(block)
    out{k} = textscan(block{k}, '%f', 'delimiter', ' ', ...
        'MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
dDot = cell2mat(out);

% Read in the vector of the acceleration DOFs
block = regexp(vtkFileNameRatesCurrentTime, 'dDDot ', 'split');
block(1) = [];
out = cell(size(block));
for k = 1:numel(block)
    out{k} = textscan(block{k}, '%f', 'delimiter', ' ', ...
        'MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
dDDot = cell2mat(out);

end