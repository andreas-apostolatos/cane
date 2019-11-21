function massMtx = computeMassMtx4FEMVMSStabNSE2D...
    (analysis,fldMsh,parameters,gaussInt)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the mass matrix correponding to a finite element formulation of 
% the 2D Navier-Stokes problem.
%
%             Input :
%          analysis : Information on the analysis
%                           .type : Analysis type
%            fldMsh : Nodes and elements of the fluid mesh
%        parameters : Flow parameters
%          gaussInt : Structure responsible for the integration
%
%            Output :
%           massMtx : The mass matrix of the Navier-Stokes problem
%
% Function layout :
%
%% Function main body

%% 0. Read input

% Initialize arrays
massMtx = 'undefined';

% This function needs to be investigated as the mass matrix seems to be
% dependent on the stabilization terms

end
