function BStrain = computeBOperatorMatrix4StrainIGAKirchhoffLoveShellLinear...
    (p,q,dR,GCovariant)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the B-operator matrix for the strain of the linear isogeometric
% Kirchhoff-Love shell formulation in the local Cartesian space.
%
%           Input : 
%             p,q : polynomial degrees
%              dR : The basis function and up to its second mixed
%                   derivatives
%      GCovariant : = [GCov1 GCov2 GCov3] the covariant base vectors
%                   derivatives of the covariant base vectors neglecting 
%                   the surface normal base vector
%
%          Output :
%         BStrain : The B-operator for the strains
%
% Function layout :
%
% 0. Read input
%
% 1. Compute the covariant metric coefficient tensor
%
% 2. Compute the contravariant base vectors
%
% 3. Compute the local Cartesian basis
%
% 4. Compute the transformation from the contravariant to the local Cartesian basis
%
% 5. Loop over all the degrees of freedom to compute the B-operator matrix in the curvilinear system
% ->
%    5i. Compute local node number k and dof direction dir
%
%   5ii. Compute the B-operator matrix for the strain in the curvilinear space
% <-
% 6. Compute the B-operator matrix in the local Cartesian basis
%
%% Function main body

%% 0. Read input

% Compute the number of DOFs
nDOFsEl = 3*(p + 1)*(q + 1);

% Initialize auxiliary array
BStrainCurvilinear = zeros(3,nDOFsEl);

%% 1. Compute the covariant metric coefficient tensor
GabCov = GCovariant'*GCovariant;

%% 2. Compute the contravariant base vectors
GContravariant = (GabCov\GCovariant')';

%% 3. Compute the local Cartesian basis
eLC = computeLocalCartesianBasis4BSplineSurface(GCovariant,GContravariant);

%% 4. Compute the transformation from the contravariant to the local Cartesian basis
TContra2LCVoigt = computeTFromContra2LocalCartesian4VoigtStrainIGAKLShell(eLC,GContravariant);

%% 5. Loop over all the degrees of freedom to compute the B-operator matrix in the curvilinear system
for r = 1:nDOFsEl
    %% 5i. Compute local node number k and dof direction dir
    k = ceil(r/3);
    dir = r - 3*(k-1);

    %% 5ii. Compute the B-operator matrix for the strain in the curvilinear space
    BStrainCurvilinear(1,r) = dR(k,2)*GCovariant(dir,1);
    BStrainCurvilinear(2,r) = dR(k,4)*GCovariant(dir,2);
    BStrainCurvilinear(3,r) = 1/2*(dR(k,2)*GCovariant(dir,2) + GCovariant(dir,1)*dR(k,4));
end

%% 6. Compute the B-operator matrix in the local Cartesian basis
BStrain = TContra2LCVoigt*BStrainCurvilinear;

end

