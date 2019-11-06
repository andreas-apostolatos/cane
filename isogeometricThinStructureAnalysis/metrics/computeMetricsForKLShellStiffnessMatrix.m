function [GCovariant,dGCovariant,GabCovariant,GContravariant,G3Tilde,CurvatureVoigt,dA,G3,eLC,TFromContraToLC] = computeMetricsForKLShellStiffnessMatrix(i,p,u,U,j,q,v,V,CP,dR,ddR)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the covariant and the contravariant basis, the metric 
% coefficients for the covariant basis, the contravariant basis, the not
% normalized orthogonal to the surface base vector, the differential 
% element area as well as the transformation matrix from the contravariant 
% basis to the local Cartesian one
%
%           Input :
%             i,j : knot span indices
%             p,q : Polynomial degrees
%             u,v : parametric coordinates on the surfaces
%             U,V : knot spans in u,v-directions
%              CP : control point coordinates and weights
%              dR : The first derivatives of the basis functions
%             ddR : The second derivatives of the basis functions
%
%          Output :
%      GCovariant : The base vectors of the 2D surface g_a = g(:,a) a = 1,2
%     dGCovariant : The derivatives of the base vectors (three columns 
%                   dg1/du dg2/dv and dg1/dv = dg2/du)
%    GabCovariant : The coefficients of the covariant metric tensor
%  GContravariant : The contravariant basis
%         G3Tilde : The not normalized normal to the surface vector
%  CurvatureVoigt : The curvatures coefficients in a Voigt notation
%              dA : The differential element area
%              G3 : The normalized normal to the surface vector (base 
%                   vector G3)
%             eLC : The local Cartesian basis, namely, 
%                   eLC1 = GCov1/norm(GCov1),
%                   eLC2 = GContra2/norm(GContra2),
%                   eLC3 = G3
% TFromContraToLC : The transformation matrix from the contravariant to the
%                   local Cartesian basis
%
% Function layout :
%
% 1. Compute the covariant base vectors g and their derivatives
%
% 2. Compute the normal to the surface vector G3Tilde (not normalized)
%
% 3. Compute the legth of G3 (= area dA)
%
% 4. Compute normal vector G3 (third base vector for both the covariant and the contravariant base)
%
% 5. Compute curvature coefficients arranged in Voigt notation
%
% 6. Compute the covariant metric tensor G_ab
%
% 7. Compute the contravariant metric tensor G^ab via (G_ab)^-1
%
% 8. Compute the local Cartesian basis
%
% 9. Compute the transformation matrix from the contravariant to local cartesian basis
% 
%% Function documentation

%% 1. Compute the covariant base vectors g and their derivatives
[GCovariant,dGCovariant] = computeBaseVectorsAndDerivatives2DGivenBasisFunctions(i,p,u,U,j,q,v,V,CP,dR,ddR);

%% 2. Compute the normal to the surface vector G3Tilde (not normalized)
G3Tilde = cross(GCovariant(:,1),GCovariant(:,2));

%% 3. Compute the legth of G3 (= area dA)
dA = norm(G3Tilde);

%% 4. Compute normal vector G3 (third base vector for both the covariant and the contravariant base)
G3 = G3Tilde/dA;

%% 5. Compute curvature coefficients arranged in Voigt notation
CurvatureVoigt = dGCovariant'*G3;

%% 6. Compute the covariant metric tensor G_ab
GabCovariant = zeros(2,2);
GabCovariant(1,1) = GCovariant(:,1)'*GCovariant(:,1);
GabCovariant(1,2) = GCovariant(:,1)'*GCovariant(:,2);
GabCovariant(2,1) = GabCovariant(1,2);
GabCovariant(2,2) = GCovariant(:,2)'*GCovariant(:,2);

%% 7. Compute the contravariant metric tensor G^ab via (G_ab)^-1
GContravariant = GabCovariant\GCovariant';
GContravariant = GContravariant';

%% 8. Compute the local Cartesian basis

% Normalize the contravariant base vectors (G_1_hat and G^2_hat will be used 
% as the local Cartesian basis since they are orthogonal to each other, yet
% both being tangent to the surface)

% As first use the covariant base vector GCov1
eLC = zeros(3,3);
normgCov1 = norm(GCovariant(:,1));
eLC(:,1) = GCovariant(:,1)/normgCov1;

% As second use the contravariant base vector GContra2
normgContra2 = norm(GContravariant(:,2));
eLC(:,2) = GContravariant(:,2)/normgContra2;

% As third base vector use the normalized normal to the surface vector G3
eLC(:,3) = G3;

%% 9. Compute the transformation matrix from the contravariant to local cartesian basis

% Compute auxiliary scalar products between the 2 basis. This is needed
% because of the re-arrangement that the tensorial matrix has to obtain so
% that it complies with the Voigt notation of the strain and stress vectors
TFromContraToLC11 = eLC(1,1)*GContravariant(1,1) + eLC(2,1)*GContravariant(2,1) + eLC(3,1)*GContravariant(3,1);
TFromContraToLC12 = eLC(1,1)*GContravariant(1,2) + eLC(2,1)*GContravariant(2,2) + eLC(3,1)*GContravariant(3,2);
TFromContraToLC21 = eLC(1,2)*GContravariant(1,1) + eLC(2,2)*GContravariant(2,1) + eLC(3,2)*GContravariant(3,1);
TFromContraToLC22 = eLC(1,2)*GContravariant(1,2) + eLC(2,2)*GContravariant(2,2) + eLC(3,2)*GContravariant(3,2);

% Compute the entries of the transformation matrix
TFromContraToLC = zeros(3,3);
TFromContraToLC(1,1) = TFromContraToLC11*TFromContraToLC11;
TFromContraToLC(1,2) = TFromContraToLC12*TFromContraToLC12;
TFromContraToLC(1,3) = 2*TFromContraToLC11*TFromContraToLC12;
TFromContraToLC(2,1) = TFromContraToLC21*TFromContraToLC21;
TFromContraToLC(2,2) = TFromContraToLC22*TFromContraToLC22;
TFromContraToLC(2,3) = 2*TFromContraToLC21*TFromContraToLC22;
TFromContraToLC(3,1) = 2*TFromContraToLC11*TFromContraToLC21;
TFromContraToLC(3,2) = 2*TFromContraToLC12*TFromContraToLC22;
TFromContraToLC(3,3) = 2*(TFromContraToLC11*TFromContraToLC22 + TFromContraToLC12*TFromContraToLC21);

end
