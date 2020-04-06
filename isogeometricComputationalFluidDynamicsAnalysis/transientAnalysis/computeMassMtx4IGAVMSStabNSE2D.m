function massMtx = computeMassMtx4IGAVMSStabNSE2D...
    (BSplinePatches, numDOFs)
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
%    BSplinePatches : Array of the B-Spline patches comprising the
%                     multipatch geometry
%           numDOFs : Number of DOFs
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
