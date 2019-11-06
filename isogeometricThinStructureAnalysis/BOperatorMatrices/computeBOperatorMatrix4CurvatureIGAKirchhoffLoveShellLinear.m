function BCurvature = computeBOperatorMatrix4CurvatureIGAKirchhoffLoveShellLinear...
    (p,q,dR,GCovariant,dGCovariant,G3Tilde)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the B-operator matrix for the change in the curvature of the 
% linear isogeometric Kirchhoff-Love shell formulation in the local Cartesian space.
%
%           Input : 
%             p,q : polynomial degrees
%              dR : The basis function and up to its second mixed
%                   derivatives
%      GCovariant : = [GCov1 GCov2 GCov3] the covariant base vectors
%     dGCovariant : = [dGCov1/dxi dGCov2/deta dGCov1/deta = dGCov2/dxi] the
%                   firts derivatives of the covariant base vectors
%                   neglecting the surface normal
%         G3Tilde : The not normalized surface normal
%
%          Output :
%      BCurvature : The B-operator for the change in curvature
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
% 5. Loop over all the degrees of freedom to compute the B-operator matrix for the curvature change in the curvilinear system
% ->
%    5i. Compute local node number k and dof direction dir
%
%   5ii. Compute the curvature constituents
%
%  5iii. Compute the B-operator matrix for curvature in the curvilinear space Bm(:,r) = -( ddR(k,:)'*G3(dir) + dGCovariant'*dn );
% <-
% 6. Compute the B-operator matrix for the curvature change in the local Cartesian basis
%
%% Function main body

%% 0. Read input

% The number of DOFs
nDOFsEl = 3*(p + 1)*(q + 1);

% Initialize auxiliary array
BCurvatureCurvilinear = zeros(3,nDOFsEl);

%% 1. Compute the covariant metric coefficient tensor
GabCov = GCovariant'*GCovariant;

%% 2. Compute the contravariant base vectors
GContravariant = (GabCov\GCovariant')';

%% 3. Compute the local Cartesian basis
eLC = computeLocalCartesianBasis4BSplineSurface(GCovariant,GContravariant);

%% 4. Compute the transformation from the contravariant to the local Cartesian basis
TContra2LCVoigt = computeTFromContra2LocalCartesian4VoigtStrainIGAKLShell(eLC,GContravariant);

%% 5. Loop over all the degrees of freedom to compute the B-operator matrix for the curvature change in the curvilinear system
for r = 1:nDOFsEl
    %% 5i. Compute local node number k and dof direction dir
    k = ceil(r/3);
    dir = r-3*(k-1);
    dg = zeros(3,2);
    dg(dir,1) = dR(k,2);
    dg(dir,2) = dR(k,4);

    %% 5ii. Compute the curvature constituents
    dg3 = cross(dg(:,1),GCovariant(:,2)) + cross(GCovariant(:,1),dg(:,2));
    g3dg3lg3_3 = G3Tilde'*dg3/norm(G3Tilde)^3;
    dn = dg3/norm(G3Tilde) - G3Tilde*g3dg3lg3_3;

    %% 5iii. Compute the B-operator matrix for curvature in the curvilinear space Bm(:,r) = -( ddR(k,:)'*G3(dir) + dGCovariant'*dn );
    BCurvatureCurvilinear(:,r) = - (dR(k,[3 6 5])'*G3Tilde(dir)/norm(G3Tilde) + dGCovariant'*dn);
end

%% 6. Compute the B-operator matrix for the curvature change in the local Cartesian basis
BCurvature = TContra2LCVoigt*BCurvatureCurvilinear;

end

