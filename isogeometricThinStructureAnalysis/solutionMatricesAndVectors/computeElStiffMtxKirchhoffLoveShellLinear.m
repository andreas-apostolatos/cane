function Ke = computeElStiffMtxKirchhoffLoveShellLinear ...
    (p, q, dR, GCovariant, dGCovariant, G3Tilde, Dm, Db)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the element stiffness matrix corresponding to the linear
% isogeometric Kirchhoff-Love shell problem.
%
%           Input : 
%             p,q : polynomial degrees
%              dR : The basis function and up to its second mixed
%                   derivatives
%      GCovariant : = [GCov1 GCov2 GCov3] the covariant base vectors
%     dGCovariant : = [dGCov1/dxi dGCov2/deta dGCov1/deta=dGCov2/dxi] the 
%                   derivatives of the covariant base vectors neglecting 
%                   the surface normal base vector
%         G3Tilde : The not normalized surface normal
%           Dm,Db : membrane and bending material matrices
%
%          Output :
%              ke : The element stiffness matrix at the Gauss Point
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
% 4. Compute the transformation matrix from the contravariant to the local Cartesian basis
%
% 5. Compute the B-operator matrices for the mebrane and bending part of the stiffnesses in the curvilinear system
%
% 6. Compute the B-operator matrices for the mebrane and bending part of the stiffnesses in the local Cartesian system
%
% 7. Compute the element stiffness matrices for the membrane and the bending part
%
% 8. Add both contributions from the membrane as well as the bending part to the element stiffness matrix
%
%% Function main body

%% 0. Read input

% The number of DOFs
numCPsEl = (p + 1)*(q + 1);
numDOFsEl = 3*numCPsEl;

% Initialize auxiliary arrays
BOperatorMatrixMembraneCurvilinear = zeros(3,numDOFsEl);
BOperatorMatrixBendingCurvilinear = zeros(3,numDOFsEl);

%% 1. Compute the covariant metric coefficient tensor
GabCov = GCovariant'*GCovariant;

%% 2. Compute the contravariant base vectors
GContravariant = (GabCov\GCovariant')';

%% 3. Compute the local Cartesian basis
eLC = computeLocalCartesianBasis4BSplineSurface ...
    (GCovariant, GContravariant);

%% 4. Compute the transformation matrix from the contravariant to the local Cartesian basis
TContra2LC = computeTFromContra2LocalCartesian4VoigtStrainIGAKLShell ...
    (eLC, GContravariant);

%% 5. Compute the B-operator matrices for the mebrane and bending part of the stiffnesses in the curvilinear system

% Loop over all the degrees of freedom
for iDOFs = 1:numDOFsEl
    % Compute local node number k and dof direction dir
    k = ceil(iDOFs/3);
    dir = iDOFs - 3*(k - 1);
    dg = zeros(3, 2);
    dg(dir, 1) = dR(k, 2);
    dg(dir, 2) = dR(k, 4);

    % Compute the B-operator matrix for the membrane part in the
    % curvilinear space
    BOperatorMatrixMembraneCurvilinear(1, iDOFs) = dR(k, 2)*GCovariant(dir, 1);
    BOperatorMatrixMembraneCurvilinear(2, iDOFs) = dR(k, 4)*GCovariant(dir, 2);
    BOperatorMatrixMembraneCurvilinear(3, iDOFs) = 1/2*(dR(k, 2)*GCovariant(dir, 2) + ...
        GCovariant(dir, 1)*dR(k, 4));

    % Compute the curvature constituents
    dg3 = cross(dg(:, 1), GCovariant(:, 2)) + cross(GCovariant(:, 1), dg(:, 2));
    g3dg3lg3_3 = G3Tilde'*dg3/norm(G3Tilde)^3;
    dn = dg3/norm(G3Tilde) - G3Tilde*g3dg3lg3_3;

    % Compute the B-operator matrix for the bending part in the curvilinear 
    % space
    BOperatorMatrixBendingCurvilinear(:, iDOFs) = ...
        -(dR(k, [3 6 5])'*G3Tilde(dir)/norm(G3Tilde) + dGCovariant'*dn);
end

%% 6. Compute the B-operator matrices for the mebrane and bending part of the stiffnesses in the local Cartesian system

% For the membrane part
BOperatorMatrixMembraneCartesian = TContra2LC*BOperatorMatrixMembraneCurvilinear;

% For the bending part
BOperatorMatrixBendingCartesian = TContra2LC*BOperatorMatrixBendingCurvilinear;

%% 7. Compute the element stiffness matrices for the membrane and the bending part

% Compute the membrane stiffness
Kem = BOperatorMatrixMembraneCartesian'*Dm*BOperatorMatrixMembraneCartesian;

% Compute the bending stiffness
Keb = BOperatorMatrixBendingCartesian'*Db*BOperatorMatrixBendingCartesian;

%% 8. Add both contributions from the membrane as well as the bending part to the element stiffness matrix
Ke = Kem + Keb;

end