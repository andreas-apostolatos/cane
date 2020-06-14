function [BStrainLC, BCurvatureLC, dBCurvatureLC1, dBCurvatureLC2] = ...
    computeBOperatorMatrices4IGAKirchhoffLoveShellLinear ...
    (p, q, dR, GCovariant, dGCovariant, ddGCovariant, G3Tilde, TContra2LC)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the B-operator matrices for the strain and bending parts of the
% stiffness as well as the B-operator matrices of the derivatives of the
% change in the curvature tensors in both parametric directions in the 
% local Cartesian basis for the isogeometric Kirchhoff-Love shell.
%
%           Input : 
%             p,q : polynomial degrees
%              dR : The basis function and up to its third mixed
%                   derivatives
%      GCovariant : = [GCov1 GCov2 GCov3] the covariant base vectors
%     dGCovariant : = [dGCov1/dxi dGCov2/deta dGCov1/deta = dGCov2/dxi] the
%                   firts derivatives of the covariant base vectors 
%    ddGCovariant : = [d^2GCov1/dxi^2 d^2GCov2/deta^2 d^2GCov1/dxi*deta 
%                     d^2GCov1/deta^2] the second derivatives of the
%                     covariant base vectors
%         G3Tilde : The not normalized surface normal
%
%          Output :
%       BStrainLC : The B-operator matrix of the strains in the local
%                   Cartesian system
%    BCurvatureLC : The B-operator matrix of the change in curvature in the
%                   local Cartesian system
%  dBCurvatureLC1 : The B-operator matrix of the derivatives of the
%                   curvature with respect to theta^1 covariant coordinate
%                   in the local Cartesian system
%  dBCurvatureLC2 : The B-operator matrix of the derivatives of the
%                   curvature with respect to theta^2 covariant coordinate
%                   in the local Cartesian system
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the degrees of freedom to compute the B-operator matrix in the curvilinear system
% ->
%    1i. Compute local node number k and dof direction dir and the variation of the derivatives of the covariant base vectors
%
%   1ii. Compute the curvature constituents
%
%  1iii. Compute variations of the base vectors and their derivatives needed for the computation of the derivatives of the curvature variation
%
%   1iv. Compute the B-operator matrix for the strain in contravariant basis
%
%    1v. Compute the B-operator matrix for curvature in contravariant basis Bm(:,r) = -( ddR(k,:)'*G3(dir) + dGCovariant'*dn );
%
%   1vi. Compute the B-operator matrix for the curvature derivatives
% <-
% 2. Compute the B-operator matrix for the strain in the local Cartesian basis
%
% 3. Compute the B-operator matrix for the curvature change in the local Cartesian basis
%
% 4. Compute the B-operator matrix for the derivatives of the the curvature with respect to theta^1 and theta^2 directions in the local Cartesian system
%
%% Function main body

%% 0. Read input

% Compute the number of DOFs
nDOFsEl = 3*(p + 1)*(q + 1);

% Initialize auxiliary arrays
BStrainContravariant = zeros(3,nDOFsEl);
BCurvatureContravariant = zeros(3,nDOFsEl);
dBCurvatureContravariant1 = zeros(3,nDOFsEl);
dBCurvatureContravariant2 = zeros(3,nDOFsEl);
dg = zeros(3,2);

% Preliminary computations
dg3_1 = cross(dGCovariant(:,1),GCovariant(:,2)) + cross(GCovariant(:,1),dGCovariant(:,3));
dg3_2 = cross(dGCovariant(:,3),GCovariant(:,2)) + cross(GCovariant(:,1),dGCovariant(:,2));
dn_1 = dg3_1/norm(G3Tilde) - G3Tilde*(G3Tilde'*dg3_1)/norm(G3Tilde)^3;
dn_2 = dg3_2/norm(G3Tilde) - G3Tilde*(G3Tilde'*dg3_2)/norm(G3Tilde)^3;

%% 1. Loop over all the degrees of freedom to compute the B-operator matrix in the curvilinear system
for r = 1:nDOFsEl
    %% 1i. Compute local node number k and dof direction dir and the variation of the derivatives of the covariant base vectors
    k = ceil(r/3);
    dir = r - 3*(k-1);
    dg(dir,1) = dR(k,2);
    dg(dir,2) = dR(k,5);
    
    %% 1ii. Compute the curvature constituents
    dg3 = cross(dg(:,1),GCovariant(:,2)) + cross(GCovariant(:,1),dg(:,2));
    g3dg3lg3_3 = G3Tilde'*dg3/norm(G3Tilde)^3;
    dn = dg3/norm(G3Tilde) - G3Tilde*g3dg3lg3_3;
    
    %% 1iii. Compute variations of the base vectors and their derivatives needed for the computation of the derivatives of the curvature variation
    
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
	
    %% 1iv. Compute the B-operator matrix for the strain in the contravariant basis
    BStrainContravariant(1,r) = dR(k,2)*GCovariant(dir,1);
    BStrainContravariant(2,r) = dR(k,5)*GCovariant(dir,2);
    BStrainContravariant(3,r) = 1/2*(dR(k,2)*GCovariant(dir,2) + GCovariant(dir,1)*dR(k,5));
    
    %% 1v. Compute the B-operator matrix for curvature in the contravariant basis Bm(:,r) = -( ddR(k,:)'*G3(dir) + dGCovariant'*dn );
    BCurvatureContravariant(:,r) = -(dR(k,[3 8 6])'*G3Tilde(dir)/norm(G3Tilde) + dGCovariant'*dn);
    
    %% 1vi. Compute the B-operator matrix for the curvature derivatives
    
    % For the derivatives w.r.t. coordinate theta^1
    dBCurvatureContravariant1(1,r) = -(dR(k,4)*G3Tilde(dir)/norm(G3Tilde) + dR(k,3)*dn_1(dir) + ddGCovariant(:,1)'*dn + dGCovariant(:,1)'*ddn_1);
    dBCurvatureContravariant1(2,r) = -(dR(k,9)*G3Tilde(dir)/norm(G3Tilde) + dR(k,8)*dn_1(dir) + ddGCovariant(:,4)'*dn + dGCovariant(:,2)'*ddn_1);
    dBCurvatureContravariant1(3,r) = -(dR(k,7)*G3Tilde(dir)/norm(G3Tilde) + dR(k,6)*dn_1(dir) + ddGCovariant(:,3)'*dn + dGCovariant(:,3)'*ddn_1);
    
    % For the derivatives w.r.t. coordinate theta^2
    dBCurvatureContravariant2(1,r) = -(dR(k,7)*G3Tilde(dir)/norm(G3Tilde) + dR(k,3)*dn_2(dir) + ddGCovariant(:,3)'*dn + dGCovariant(:,1)'*ddn_2);
    dBCurvatureContravariant2(2,r) = -(dR(k,10)*G3Tilde(dir)/norm(G3Tilde) + dR(k,8)*dn_2(dir) + ddGCovariant(:,2)'*dn + dGCovariant(:,2)'*ddn_2);
    dBCurvatureContravariant2(3,r) = -(dR(k,9)*G3Tilde(dir)/norm(G3Tilde) + dR(k,6)*dn_2(dir) + ddGCovariant(:,4)'*dn + dGCovariant(:,3)'*ddn_2);
end

%% 2. Compute the B-operator matrix for the strain in the local Cartesian basis
BStrainLC = TContra2LC*BStrainContravariant;

%% 3. Compute the B-operator matrix for the curvature change in the local Cartesian basis
BCurvatureLC = TContra2LC*BCurvatureContravariant;

%% 4. Compute the B-operator matrix for the derivatives of the the curvature with respect to theta^1 and theta^2 directions in the local Cartesian system

% Derivatives w.r.t. theta^1 direction
dBCurvatureLC1 = TContra2LC*dBCurvatureContravariant1;

% Derivatives w.r.t. theta^2 direction
dBCurvatureLC2 = TContra2LC*dBCurvatureContravariant2;

end
