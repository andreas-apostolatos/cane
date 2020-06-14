function massMtx = computeMassMtx4FEMVMSStabNSE2D...
    (propAnalysis, fldMsh, parameters, propGaussInt)
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
%      propAnalysis : Structure containing general information on the
%                     analysis,
%                           .type : Analysis type
%            fldMsh : Nodes and elements of the fluid mesh
%        parameters : Flow parameters
%      propGaussInt : Structure containing information on the numerical
%                     quadrature
%
%            Output :
%           massMtx : The mass matrix of the Navier-Stokes problem
%
%% Function main body

% Initialize arrays
massMtx = 'undefined';

% This function needs to be investigated as the mass matrix is dependent on 
% the stabilization terms

end