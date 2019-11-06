function strainVoigt = computePostprocVoigtStrainIGAKirchhoffLoveShellNLinear...
    (GCovariant,gCovariant)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the strain vector in a Voigt notation for the postprocessing of
% the nonlinear analysis performed over the isogeometric Kirchhoff-Love
% shell.
%
%           Input : 
%      GCovariant : = [GCov1 GCov2 GCov3] the covariant base vectors and
%                   their derivatives neglecting the derivatives  of the 
%                   surface normal base vector for the reference
%                   configuration
%      gCovariant :  = [gCov1 gCov2 gCov3] the covariant base vectors and
%                   their derivatives neglecting the derivatives  of the 
%                   surface normal base vector for the current
%                   configuration
%
%          Output :
%     strainVoigt : The nonlinear strain vector 
%                   [epsilon_11 epsilon_22 epsilon_12]
%
% Function layout :
%
% 0. Read input
%
% 1. Compute metrics for the reference and the current configuration
%
% 2. Compute the strain in a Voigt notation in the curvilinear space
%
% 3. Compute the strain in a Voigt notation in the local Cartesian space
%
%% Function main body

%% 0. Read input

% Initialize arrays
strainVoigtCurvilinear = zeros(3,1);

%% 1. Compute metrics for the reference and the current configuration

% For the reference configuration :
% _________________________________

% Compute the covariant metric coefficient tensor
GabCov = GCovariant'*GCovariant;

% Compute the contravariant basis
GContravariant = (GabCov\GCovariant')';

% Compute the local Cartesian basis
eLC = computeLocalCartesianBasis4BSplineSurface(GCovariant,GContravariant);

% Compute the transformation from the contravariant to the local Cartesian 
% basis for a Voigt-type second order tensor
TContra2LC = computeTFromContra2LocalCartesian4VoigtStrainIGAKLShell(eLC,GContravariant);

% For the current configuration :
% _________________________________

% Compute the covariant metric coefficient tensor
gabCov = gCovariant'*gCovariant;

%% 2. Compute the strain in a Voigt notation in the curvilinear space
strainVoigtCurvilinear(1:2,1) = 1/2*(gabCov(1:2,1)-GabCov(1:2,1));
strainVoigtCurvilinear(3,1) = 1/2*(gabCov(1,2)-GabCov(1,2));

%% 3. Compute the strain in a Voigt notation in the local Cartesian space
strainVoigt = TContra2LC*strainVoigtCurvilinear;

end
