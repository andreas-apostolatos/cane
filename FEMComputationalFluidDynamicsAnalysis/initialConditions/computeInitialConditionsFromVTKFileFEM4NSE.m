function [up, upDot, upDDot, numTimeStep] = ...
    computeInitialConditionsFromVTKFileFEM4NSE...
    (analysis, fldMsh, DOF4Output, parameters, fldDynamics, ...
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
% Returns the initial velocity and pressure fields as well as their rates,
% for a CFD simulation which has been stopped, by reading the results from
% the generated output files defined in VTKResultFile.
%
%             Input :
%          analysis : On the analysis :
%                       .type : Analysis type
%            fldMsh : Nodes and elements of the fluid mesh
%        DOF4Output : Array containing the arrangment of the DOFs for 
%                     printing them out
%        parameters : Flow parameters (null variable for this function)
%       fldDynamics : Transient analysis parameters for the CFD simulation
%     VTKResultFile : The name of the result file in the output folder
%                     where to get the initial conditions for the transient 
%                     simulation
%          caseName : The name of the case
%        pathToFile : Path to the file where to get the initial conditions 
%                     from
%
%            Output :
%                up : The vector of DoFs containing the initial conditions 
%                     for the velocity and pressure field
%             upDot : The vector of DoFs containing the initial conditions 
%                     for the acceleration and the pressure rate field
%            upDDot : Dummny variable needed only for computational 
%                     structural dynamics
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
% 3. Get the velocity and the pressure field of the current time step
%
% 4. Get the rates of the velocities and the pressure fields of the current time step
%
%% Function main body

%% 0. Read input

% Number of nodes
noNodes = length(fldMsh.nodes(:,1));

% Number of degees of freedom
noDOFs = 3*noNodes;

% Get the time step when the simulation was stopped
[~,vtkFileName,~] = fileparts(strcat(pathToFile,VTKResultFile,'.vtk'));

% Initialize output arrays
up = zeros(noDOFs,1);
upDDot = 'undefined';

%% 1. Extract the time step number from the current and the previous time step
counter = 1;
for i=length(vtkFileName):-1:1
    if ~strcmp(vtkFileName(i),'_')
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

%% 3. Get the velocity and the pressure field of the current time step

% Read the vtk output file for the current time step
vtkFileNameCurrentTimeStep = fileread(strcat(pathToFile,caseName,'/',caseName,VTKResultFile,'.vtk'));

% Read in the velocity field
block = regexp(vtkFileNameCurrentTimeStep,'VECTORS velocity double','split');
block(1) = [];
out = cell(size(block));
for k = 1:numel(block)
    out{k} = textscan(block{k},'%f %f %f','delimiter',' ','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
velocities = cell2mat(out);

% Read in the pressure field
block = regexp(vtkFileNameCurrentTimeStep,'SCALARS pressure double\nLOOKUP_TABLE default','split');
block(1) = [];
out = cell(size(block));
for k = 1:numel(block)
    out{k} = textscan(block{k},'%f','delimiter',' ','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
pressures = cell2mat(out);

% Assemble the vector of the initial conditions for the unknown vector
% field up
for i=1:2
    up(DOF4Output(i,:)) = velocities(:,i);
end
up(DOF4Output(3,:)) = pressures;

%% 4. Get the rates of the velocities and the pressure fields of the current time step

% Get the file name of the VTK result file for the rates of the current
% time step
VTKResultFileRatesCurrentTimeStep = strcat('_upRates_',timeStepNoString);

% Read the file
vtkFileNameRatesCurrentTime = ...
    fileread(strcat(pathToFile,caseName,'/upRates/',caseName,VTKResultFileRatesCurrentTimeStep,'.txt'));

% Read in the vector of the rates
block = regexp(vtkFileNameRatesCurrentTime,'upRates ','split');
block(1) = [];
out = cell(size(block));
for k=1:numel(block)
    out{k} = textscan(block{k},'%f','delimiter',' ','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
upDot = cell2mat(out);

end
