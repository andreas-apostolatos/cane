function BShearForce = computeBOperatorMatrix4ShearForcesIGAKirchhoffLoveShellLinear(p,q,dR,GCovariant,dGCovariant,ddGCovariant,G3Tilde,Db)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the B-operator matrix for the shear forces corresponding to the
% linear isogeometric Kirchhoff-Love shell formulation in the local 
% Cartesian space. Within the Kirchhoff-Love shell theory the shear forces
% are simply the derivatives of the moments out of the equilibrium
% equations in the curvilinear system.
%
%           Input : 
%             p,q : polynomial degrees
%              dR : The basis function and up to its third mixed
%                   derivatives
%      GCovariant : = [GCov1 GCov2 GCov3] the covariant base vectors
%                   derivatives of the covariant base vectors neglecting 
%                   the surface normal base vector
%     dGCovariant : = [dGCov1/dxi dGCov2/deta dGCov1/deta = dGCov2/dxi] the
%                   first derivatives of the covariant base vectors
%    ddGCovariant : = [d^2GCov1/dxi^2 d^2GCov2/deta^2 d^2GCov1/dxi*deta 
%                     d^2GCov1/deta^2] the second derivatives of the
%                     covariant base vectors
%         G3Tilde : The not normalized surface normal
%              Db : The material matrix of the bending part of the
%                   stiffness
%
%          Output :
%     BShearForce : The B-operator matrix for the shear forces 2xnDOFsEl
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
% 5. Loop over all the degrees of freedom to compute the B-operator matrix for the derivatives of the curvature with respect to theta^1 and theta^2 directions in the curvilinear system
% ->
%    5i. Compute local node number k and dof direction dir
%
%   5ii. Compute variations of the base vectors and their derivatives needed for the computation of the curvature variation
%
%  5iii. Compute variations of the base vectors and their derivatives needed for the computation of the derivatives of the curvature variation
%
%   5iv. Compute the B-operator matrix for the curvature derivatives
% <-
% 6. Compute the B-operator matrix for the derivatives of the the curvature with respect to theta^1 and theta^2 directions in the local Cartesian system
%
% 7. Compute the B-operator matrices for the shear forces
%
%% Function main body

%% 0. Read input

% The number of DOFs
nDOFsEl = 3*(p + 1)*(q + 1);

% Initialize auxiliary arrays
dg = zeros(3,2);
dg3 = zeros(3,1);
dn = zeros(3,1);
dBCurvatureCurvilinear1 = zeros(3,nDOFsEl);
dBCurvatureCurvilinear2 = zeros(3,nDOFsEl);

% Preliminary computations
dg3_1 = cross(dGCovariant(:,1),GCovariant(:,2)) + cross(GCovariant(:,1),dGCovariant(:,3));
dg3_2 = cross(dGCovariant(:,3),GCovariant(:,2)) + cross(GCovariant(:,1),dGCovariant(:,2));
dn_1 = dg3_1/norm(G3Tilde) - G3Tilde*(G3Tilde'*dg3_1)/norm(G3Tilde)^3;
dn_2 = dg3_2/norm(G3Tilde) - G3Tilde*(G3Tilde'*dg3_2)/norm(G3Tilde)^3;


%% 1. Compute the covariant metric coefficient tensor
GabCov = GCovariant'*GCovariant;

%% 2. Compute the contravariant base vectors
GContravariant = (GabCov\GCovariant')';

%% 3. Compute the local Cartesian basis
eLC = computeLocalCartesianBasis4BSplineSurface(GCovariant,GContravariant);

%% 4. Compute the transformation from the contravariant to the local Cartesian basis
TContra2LCVoigt = computeTFromContra2LocalCartesian4VoigtStrainIGAKLShell(eLC,GContravariant);

%% 5. Loop over all the degrees of freedom to compute the B-operator matrix for the derivatives of the curvature with respect to theta^1 and theta^2 directions in the curvilinear system
for r = 1:nDOFsEl
    %% 5i. Compute local node number k and dof direction dir
    k = ceil(r/3);
    dir = r-3*(k-1);
    dg(dir,1) = dR(k,2);
    dg(dir,2) = dR(k,5);

    %% 5ii. Compute variations of the base vectors and their derivatives needed for the computation of the curvature variation
    dg3(1) = dg(2,1)*GCovariant(3,2)-dg(3,1)*GCovariant(2,2) + GCovariant(2,1)*dg(3,2)-GCovariant(3,1)*dg(2,2);
    dg3(2) = dg(3,1)*GCovariant(1,2)-dg(1,1)*GCovariant(3,2) + GCovariant(3,1)*dg(1,2)-GCovariant(1,1)*dg(3,2);
    dg3(3) = dg(1,1)*GCovariant(2,2)-dg(2,1)*GCovariant(1,2) + GCovariant(1,1)*dg(2,2)-GCovariant(2,1)*dg(1,2);

    g3dg3lg3_3 = G3Tilde'*dg3/norm(G3Tilde)^3;

    dn(1) = dg3(1)/norm(G3Tilde) - G3Tilde(1)*g3dg3lg3_3;
    dn(2) = dg3(2)/norm(G3Tilde) - G3Tilde(2)*g3dg3lg3_3;
    dn(3) = dg3(3)/norm(G3Tilde) - G3Tilde(3)*g3dg3lg3_3;
    
    %% 5iii. Compute variations of the base vectors and their derivatives needed for the computation of the derivatives of the curvature variation
    
    % For the derivatives w.r.t. coordinate theta^1
    ddg3_1 = cross(ddGCovariant(:,1),GCovariant(:,2)) + cross(dg(:,1),dGCovariant(:,3)) ...
         + cross(dGCovariant(:,1),dg(:,2)) + cross(GCovariant(:,1),ddGCovariant(:,3));                
    ddn_1 = ddg3_1/norm(G3Tilde) - (dg3*(G3Tilde'*dg3_1) + dg3_1*(G3Tilde'*dg3) ...
         + G3Tilde*(ddg3_1'*G3Tilde + dg3'*dg3_1))/norm(G3Tilde)^3 ...
         + 3*G3Tilde*(G3Tilde'*dg3)*(G3Tilde'*dg3_1)/norm(G3Tilde)^5;  
   
    % For the derivatives w.r.t. coordinate theta^2
	ddg3_2 = cross(ddGCovariant(:,3),GCovariant(:,2)) + cross(dg(:,1),dGCovariant(:,2)) ...
         + cross(dGCovariant(:,3),dg(:,2)) + cross(GCovariant(:,1),ddGCovariant(:,2));                
    ddn_2 = ddg3_2/norm(G3Tilde) - (dg3*(G3Tilde'*dg3_2) + dg3_2*(G3Tilde'*dg3) ...
         + G3Tilde*(ddg3_2'*G3Tilde + dg3'*dg3_2))/norm(G3Tilde)^3 ...
         + 3*G3Tilde*(G3Tilde'*dg3)*(G3Tilde'*dg3_2)/norm(G3Tilde)^5;  
    
    %% 5iv. Compute the B-operator matrix for the curvature derivatives
    
    % For the derivatives w.r.t. coordinate theta^1
    dBCurvatureCurvilinear1(1,r) = -(dR(k,4)*G3Tilde(dir)/norm(G3Tilde) + dR(k,3)*dn_1(dir) + ddGCovariant(:,1)'*dn + dGCovariant(:,1)'*ddn_1);
    dBCurvatureCurvilinear1(2,r) = -(dR(k,9)*G3Tilde(dir)/norm(G3Tilde) + dR(k,8)*dn_1(dir) + ddGCovariant(:,4)'*dn + dGCovariant(:,2)'*ddn_1);
    dBCurvatureCurvilinear1(3,r) = -(dR(k,7)*G3Tilde(dir)/norm(G3Tilde) + dR(k,6)*dn_1(dir) + ddGCovariant(:,3)'*dn + dGCovariant(:,3)'*ddn_1);
    
    % For the derivatives w.r.t. coordinate theta^2
    dBCurvatureCurvilinear2(1,r) = -(dR(k,7)*G3Tilde(dir)/norm(G3Tilde) + dR(k,3)*dn_2(dir) + ddGCovariant(:,3)'*dn + dGCovariant(:,1)'*ddn_2);
    dBCurvatureCurvilinear2(2,r) = -(dR(k,10)*G3Tilde(dir)/norm(G3Tilde) + dR(k,8)*dn_2(dir) + ddGCovariant(:,2)'*dn + dGCovariant(:,2)'*ddn_2);
    dBCurvatureCurvilinear2(3,r) = -(dR(k,9)*G3Tilde(dir)/norm(G3Tilde) + dR(k,6)*dn_2(dir) + ddGCovariant(:,4)'*dn + dGCovariant(:,3)'*ddn_2);
end

%% 6. Compute the B-operator matrix for the derivatives of the the curvature with respect to theta^1 and theta^2 directions in the local Cartesian system

% Derivatives w.r.t. theta^1 direction
dBCurvatureLocalCartesian1 = TContra2LCVoigt*dBCurvatureCurvilinear1;

% Derivatives w.r.t. theta^2 direction
dBCurvatureLocalCartesian2 = TContra2LCVoigt*dBCurvatureCurvilinear2;

%% 7. Compute the B-operator matrices for the shear forces

% Shear force q^1
dShearForce1 = Db*dBCurvatureLocalCartesian1;

% Shear force q^2
dShearForce2 = Db*dBCurvatureLocalCartesian2;

% The complete B-operator matrix for the shear forces
BShearForce = [dShearForce1(1,:)/norm(GCovariant(:,1)) + dShearForce2(3,:)/norm(GCovariant(:,2))
               dShearForce1(3,:)/norm(GCovariant(:,1)) + dShearForce2(2,:)/norm(GCovariant(:,2))];

end
