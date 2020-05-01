function TContra2LC4VoigtStrain = ...
    computeTFromContra2LocalCartesian4VoigtStrainIGAKLShell ...
    (eLC, GContravariant)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the transformation matrix corresponding to second order symmetric 
% tensors in a Voigt notation from the contravariant basis to the local
% Cartesian basis [eHat1 eHat2] on a B-Spline surface. This functions has
% been implemented for the isogeometric Kirchhoff-Love shell. This
% transformation also doubles the shear entry of the Voigt type vector in
% order to comply with the internal virtual work double contraction namely:
%
%                                 epsilon:n = 
%                                             |epsilon^LC_11  |' |n^LC_11|
% epsilon^LC_{alpha,beta} n^LC_{alpha,beta} = |epsilon^LC_22  |  |n^LC_22|
%                                             |2*epsilon^LC_12|  |n^LC_12|
%
% This transformation is used when given the strain vector in the
% contravariant basis as |epsilon_11 epsilon_22 epsilon_12|' multiplication 
% with TContra2LC4VoigtStrain returns 
% |epsilon^LC_11 epsilon^LC_22 2*epsilon^LC_12|'
%
%                  Input :
%                    eLC : = [eHat1 eHat2] the local Cartesian basis
%             GCovariant : = [GContra1 GContra2 GContra3] the contravariant base 
%                          vectors of the configuration
%
%                 Output :
% TContra2LC4VoigtStrain : The transformation matrix from the contravariant 
%                          to the local Cartesian basis
%
% Function layout :
%
% 0. Read input
%
% 1. Compute the transformation matrix from the contravariant to local cartesian basis
%   
%% Function main body

%% 0. Read input

% Initialize output array
TContra2LC4VoigtStrain = zeros(3,3);

% Factor that doubles the shear entry of the Voigt type 2nd order tensor
doubleFactor = 2;

%% 1. Compute the transformation matrix from the contravariant to local cartesian basis 
TContra2LC4VoigtStrain(1,1) = (GContravariant(:,1)'*eLC(:,1))^2;
TContra2LC4VoigtStrain(1,2) = (GContravariant(:,2)'*eLC(:,1))^2;
TContra2LC4VoigtStrain(1,3) = 2*(GContravariant(:,1)'*eLC(:,1))*(GContravariant(:,2)'*eLC(:,1));
TContra2LC4VoigtStrain(2,1) = (GContravariant(:,1)'*eLC(:,2))^2;
TContra2LC4VoigtStrain(2,2) = (GContravariant(:,2)'*eLC(:,2))^2;
TContra2LC4VoigtStrain(2,3) = 2*(GContravariant(:,1)'*eLC(:,2))*(GContravariant(:,2)'*eLC(:,2));
TContra2LC4VoigtStrain(3,1) = doubleFactor*((GContravariant(:,1)'*eLC(:,1))*(GContravariant(:,1)'*eLC(:,2)));
TContra2LC4VoigtStrain(3,2) = doubleFactor*((GContravariant(:,2)'*eLC(:,1))*(GContravariant(:,2)'*eLC(:,2)));
TContra2LC4VoigtStrain(3,3) = doubleFactor*((GContravariant(:,1)'*eLC(:,1))*(GContravariant(:,2)'*eLC(:,2)) + ...
    (GContravariant(:,2)'*eLC(:,1))*(GContravariant(:,1)'*eLC(:,2)));

end

