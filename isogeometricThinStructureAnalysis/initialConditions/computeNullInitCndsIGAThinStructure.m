function [dHat,dHatDot,dHatDDot,noTimeStep,propTransientAnalysis] = ...
    computeNullInitCndsIGAThinStructure...
    (BSplinePatches,noDOFs,propTransientAnalysis,caseName,pathToOutput,...
    tab,outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns zero initial conditions for a multipatch isogeometric structure
% for the primary field, its velocity and its acceleration.
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
%              caseName : Dummy variable for this function
%          pathToOutput : Dummy variable for this function
%                   tab : Dummy variable for this function
%                outMsg : Dummy variable for this function
%
%                Output :
%                  dHat : Initialization of the primary field
%               dHatDot : Initialization of the velocity of the primary
%                         field
%              dHatDDot : Initialization of the acceleration of the primary
%                         field
%            noTimeStep : The number of the starting time step
% propTransientAnalysis : The updated transient analysis properties
%                         structure
%
% Function layout :
%
% 1. Initialize the output arrays
%
% 2. Get the number of the current time step using the start time and the time step of the simulation
%
%% Function main body

%% 1. Initialize the output arrays
dHat = zeros(noDOFs,1);
dHatDot = zeros(noDOFs,1);
dHatDDot = zeros(noDOFs,1);

%% 2. Get the number of the current time step using the start time and the time step of the simulation
noTimeStep = 0;

end
