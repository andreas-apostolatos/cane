function kappaVoigt = computePostprocChangeInCurvatureIGAKirchhoffLoveShellNLinear...
    (GCovariant,G3Tilde,dGCovariant,g3Tilde,dgCovariant)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the change in the curvature tensor corresponding to the nonlinear
% analysis of the Kirchhoff-Love shell in Voigt notation.
%
%           Input :
%      GCovariant : = [GCov1 GCov2 GCov3] the covariant base vectors of the
%                   reference configuration
%         G3Tilde : The not normalized surface normal of the reference
%                   configuration
%     dGCovariant : = [dGCov1/dxi dGCov2/deta dGCov1/deta = dGCov2/dxi] the
%                   first derivatives of the covariant base vectors
%                   neglecting the surface normal of the reference
%                   configuration
%         g3Tilde : The not normalized surface normal of the current
%                   configuration
%     dgCovariant : = [dgCov1/dxi dgCov2/deta dgCov1/deta = dgCov2/dxi] the
%                   firts derivatives of the covariant base vectors
%                   neglecting the surface normal of the current
%                   configuration
%
%          Output :
%      kappaVoigt : The nonlinear curvature tensor in a Voigt notation
%                   namely kappaVoigt = [kappa_11 kappa_22 kappa_12] 
%
% Function layout :
%
% 1. Compute the necessary metrics
%
% 2. Compute the change in the curvature tensor in the curvilinear system
%
% 3. Compute the change in the curvature tensor in the local Cartesian system
%
%% Function main body

%% 1. Compute the necessary metrics

% For the reference configuration :
% _________________________________

% Compute the covariant metric coefficients
GabCov = GCovariant'*GCovariant;

% Compute the contravariant bese vectors
GContravariant = (GabCov\GCovariant')';


% Compute the local Cartesian basis
eLC = computeLocalCartesianBasis4BSplineSurface(GCovariant,GContravariant);


% Compute the transformation from the contravariant to the local Cartesian 
% basis for the Voigt type 2nd order tensor
TContra2LCVoigt = computeTFromContra2LocalCartesian4VoigtStrainIGAKLShell(eLC,GContravariant);

% Compute the curvature tensor
BV = dGCovariant'*(G3Tilde/norm(G3Tilde));

% For the current configuration :
% _______________________________

% Compute the curvature tensor
bV = dgCovariant'*(g3Tilde/norm(g3Tilde));

%% 2. Compute the change in the curvature tensor in the curvilinear system
kappaVoigtCurvilinear = -(bV - BV);

%% 3. Compute the change in the curvature tensor in the local Cartesian system
kappaVoigt = TContra2LCVoigt*kappaVoigtCurvilinear;

end
