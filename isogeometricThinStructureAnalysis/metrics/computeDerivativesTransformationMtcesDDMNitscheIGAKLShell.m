function [A11, A22, A12, dA, dTLC2ContraVoigtStraindxi, ...
    dTLC2ContraVoigtStraindeta, Gamma] = ...
    computeDerivativesTransformationMtcesDDMNitscheIGAKLShell ...
    (dA1Covariant, dA2Covariant, AContravariant, eLC, TLC2ContraVoigtStrain)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the vectors A^{11}, A^{12}, A^{22} (see equations (4a)-(4c) in 
% Nitsche - Computation Vol. 7), the matrix dA, namely:
%
%               |A^{11},1 + A^{12},2|
%           A = |A^{12},1 + A^{22},2|
%
% (see equation (2) in Nitsche -Computation Vol.7) as well as the
% derivatives of the transformation matrix TLC2ContraVoigtStrain, namely
% TLC2ContraVoigtStrain,1 and TLC2ContraVoigtStrain,2 which transforms the
% Voigt strain from the local Cartesian system to the contravariant basis.
% The application of the output arrays is in the linear Nitsche formulation
% for the multipatch Kirchhoff-Love shell coupling and especially for the
% computation of the B-operator matrix of the shear force vector.
%
%                      Input :
%               dA1Covariant : The covariant base vector in the xi-
%                              direction and its derivatives = [A1 dA1/dxi
%                              dA1/deta(=dA2/dxi)]
%               dA2Covariant : The covariant base vector in the eta-
%                              direction and its derivatives
%                              = [A2 dA2/deta]
%             AContravariant : The contravariant base vectors
%                              = [A1Contravariant A2Contravariant]
%                        eLC : The local Cartesian base vectors
%                              = [e1LC e2LC]
%      TLC2ContraVoigtStrain : The transformation matrix from the local
%                              Cartesian to the contravariant basis for the
%                              strain measure, namely:
%                 |(T^1_1)^2     (T^2_1)^2     2*T^1_1*T^2_1              |
%       TStrain = |(T^1_2)^2     (T^2_2)^2     2*T^1_2*T^2_2              |
%                 |2*T^1_1*T^1_2 2*T^2_1*T^2_2 2*(T^1_1*T^2_2+T^2_1*T^1_2)|
%
%                     Output :
%                        A11 : See equation (4a) Nitsche - Computation Vol. 
%                              7)
%                        A22 : See equation (4b) Nitsche - Computation Vol. 
%                              7)
%                        A12 : See equation (4c) Nitsche - Computation Vol. 
%                              7)
%                         dA : The matrix defined in equation (2) Nitsche -
%                              Computation Vol. 7, namely:
%                                         |A^{11},1 + A^{12},2|
%                                    dA = |A^{12},1 + A^{22},2|
%  dTLC2ContraVoigtStraindxi : The derivative of the TLC2ContraVoigtStrain 
%                              (equation (3) Nitsche - Computation Vol. 7)
%                              in the xi-parametric direction
% dTLC2ContraVoigtStraindeta : The derivative of the TLC2ContraVoigtStrain 
%                              (equation (3) Nitsche - Computation Vol. 7)
%                              in the eta-parametric direction
%                      Gamma : The matrix containing the permutation of
%                              the Christoffel symbols (see equation (12) 
%                              in Nitsche - Computation Vol. 9)
%
% Function layout :
%
% 1. Compute the derivatives of the contravariant base vectors with respect to the parametric coordinates
%
% 2. Compute the derivatives of the local Cartesian base vectors with respect to the parametric directions
%
% 3. Form the vectors A^{11}, A^{12}, A^{22} (see equations (4a)-(4c) in Nitsche - Computation Vol. 7)
%
% 4. Compute the derivatives of all entries to the transformation tensor T^{alpha}_{beta},gamma = (A^{alpha}*eLC_{beta}),gamma
%
% 5. Form the matrix dA (see equation (2) in Nitsche - Computation Vol. 7)
%
% 6. Form the derivatives of the the transformation matrix from the contravariant basis to the global Cartesian (see equation (13) in Nitsche - Computation Vol. 7)
%
% 7. Compute the Christoffel symbols (see equations (13a)-(13h) in Nitsche - Computation Vol. 9)
%
% 8. Form the Gamma matrix containing the Christoffel symbols (see equation (12) Nitsche - Computation Vol. 9)
%
%% Function main body

%% 1. Compute the derivatives of the contravariant base vectors with respect to the parametric coordinates

% Compute covariant metric determinant A (see equation (2.2) in Nitsche - 
% Computation Vol. 1) (VERIFIED)
A =  (dA1Covariant(:,1)'*dA1Covariant(:,1))*(dA2Covariant(:,1)'*dA2Covariant(:,1)) ...
    - (dA1Covariant(:,1)'*dA2Covariant(:,1))^2;

% Compute the derivatives of the metric determinant A (see equations (2.13)
% and (2.14) in Nitsche - Computation Vol. 1)

% dA/dxi (see equations (2.13)) (VERIFIED)
dAdxi = 2*(dA1Covariant(:,1)'*dA1Covariant(:,2))*(dA2Covariant(:,1)'*dA2Covariant(:,1)) + ...
    2*(dA1Covariant(:,1)'*dA1Covariant(:,1))*(dA2Covariant(:,1)'*dA1Covariant(:,3)) - ...
    2*(dA1Covariant(:,1)'*dA2Covariant(:,1))*(dA2Covariant(:,1)'*dA1Covariant(:,2) + ...
    dA1Covariant(:,1)'*dA1Covariant(:,3));

% dA/deta (see equations (2.14)) (VERIFIED)
dAdeta = 2*(dA1Covariant(:,1)'*dA1Covariant(:,3))*(dA2Covariant(:,1)'*dA2Covariant(:,1)) + ...
    2*(dA1Covariant(:,1)'*dA1Covariant(:,1))*(dA2Covariant(:,1)'*dA2Covariant(:,2)) - ...
    2*(dA1Covariant(:,1)'*dA2Covariant(:,1))*(dA2Covariant(:,1)'*dA1Covariant(:,3) + ...
    dA1Covariant(:,1)'*dA2Covariant(:,2));

% Compute the derivatives of the contravariant base vector A^1 with respect
% to both parametric directions (see equations (2.9), (2.10) in Nitsche - 
% Computation Vol. 1)

% dA^1/dxi (equations (2.9)) (VERIFIED)
dA1Contravariantdxi = - (1/A^2)*dAdxi*((dA2Covariant(:,1)'*dA2Covariant(:,1))*...
    dA1Covariant(:,1) - (dA1Covariant(:,1)'*dA2Covariant(:,1))*dA2Covariant(:,1)) + ...
    (1/A)*(2*(dA2Covariant(:,1)'*dA1Covariant(:,3))*dA1Covariant(:,1) + ...
    (dA2Covariant(:,1)'*dA2Covariant(:,1))*dA1Covariant(:,2) - ...
    (dA2Covariant(:,1)'*dA1Covariant(:,2) + dA1Covariant(:,1)'*dA1Covariant(:,3))*...
    dA2Covariant(:,1) - (dA1Covariant(:,1)'*dA2Covariant(:,1))*dA1Covariant(:,3));

% dA^1/deta (equations (2.10)) (VERIFIED)
dA1Contravariantdeta = - (1/A^2)*dAdeta*((dA2Covariant(:,1)'*dA2Covariant(:,1))*...
    dA1Covariant(:,1) - (dA1Covariant(:,1)'*dA2Covariant(:,1))*dA2Covariant(:,1)) + ...
    (1/A)*(2*(dA2Covariant(:,1)'*dA2Covariant(:,2))*dA1Covariant(:,1) + ...
    (dA2Covariant(:,1)'*dA2Covariant(:,1))*dA1Covariant(:,3) - ...
    (dA2Covariant(:,1)'*dA1Covariant(:,3) + dA1Covariant(:,1)'*dA2Covariant(:,2))*...
    dA2Covariant(:,1) - (dA1Covariant(:,1)'*dA2Covariant(:,1))*dA2Covariant(:,2));

% Compute the derivatives of the contravariant base vector A^1 with respect
% to both parametric directions (see equations (2.13), (2.14) in Nitsche - 
% Computation Vol. 1)

% dA^2/dxi (equation (2.13)) (VERIFIED)
dA2Contravariantdxi = - (1/A^2)*dAdxi*(-(dA1Covariant(:,1)'*dA2Covariant(:,1))*...
    dA1Covariant(:,1) + (dA1Covariant(:,1)'*dA1Covariant(:,1))*dA2Covariant(:,1)) + ...
    (1/A)*(-(dA2Covariant(:,1)'*dA1Covariant(:,2) + dA1Covariant(:,1)'*dA1Covariant(:,3))*...
    dA1Covariant(:,1) - (dA1Covariant(:,1)'*dA2Covariant(:,1))*dA1Covariant(:,2) + ...
    2*(dA1Covariant(:,1)'*dA1Covariant(:,2))*dA2Covariant(:,1) + ...
    (dA1Covariant(:,1)'*dA1Covariant(:,1))*dA1Covariant(:,3));

% dA^2/deta (equation (2.14)) (VERIFIED)
dA2Contravariantdeta = - (1/A^2)*dAdeta*(-(dA1Covariant(:,1)'*dA2Covariant(:,1))*...
    dA1Covariant(:,1) + (dA1Covariant(:,1)'*dA1Covariant(:,1))*dA2Covariant(:,1)) + ...
    (1/A)*(-(dA2Covariant(:,1)'*dA1Covariant(:,3) + dA1Covariant(:,1)'*dA2Covariant(:,2))*...
    dA1Covariant(:,1) - (dA1Covariant(:,1)'*dA2Covariant(:,1))*dA1Covariant(:,3) + ...
    2*(dA1Covariant(:,1)'*dA1Covariant(:,3))*dA2Covariant(:,1) + ...
    (dA1Covariant(:,1)'*dA1Covariant(:,1))*dA2Covariant(:,2));

%% 2. Compute the derivatives of the local Cartesian base vectors with respect to the parametric directions

% Compute the derivatives of the local Cartesian base vector eLC1 with
% respect to both parametric directions (see equations (5a), (5b) in 
% Nitsche - Computation Vol. 2)

% deLC^1/dxi (equation (5a))
de1LCdxi = - (1/norm(dA1Covariant(:,1))^3)*(dA1Covariant(:,1)'*dA1Covariant(:,2))*...
    dA1Covariant(:,1) + (1/norm(dA1Covariant(:,1)))*dA1Covariant(:,2);

% deLC^1/deta (equation (5b))
de1LCdeta = - (1/norm(dA1Covariant(:,1))^3)*(dA1Covariant(:,1)'*dA1Covariant(:,3))*...
    dA1Covariant(:,1) + (1/norm(dA1Covariant(:,1)))*dA1Covariant(:,3);

% Compute the derivatives of the local Cartesian base vector eLC2 with
% respect to both parametric directions (see equations (5c), (5d) in 
% Nitsche - Computation Vol. 2)

% deLC^1/dxi (equation (5c))
de2LCdxi = - (1/norm(AContravariant(:,2))^3)*(AContravariant(:,2)'*dA2Contravariantdxi)*...
    AContravariant(:,2) + (1/norm(AContravariant(:,2)))*dA2Contravariantdxi;

% deLC^1/deta (equation (5d))
de2LCdeta = - (1/norm(AContravariant(:,2))^3)*(AContravariant(:,2)'*dA2Contravariantdeta)*...
    AContravariant(:,2) + (1/norm(AContravariant(:,2)))*dA2Contravariantdeta;

%% 3. Form the vectors A^{11}, A^{12}, A^{22} (see equations (4a)-(4c) in Nitsche - Computation Vol. 7)

% Vector A^{11} (equation (4a))
A11 = [TLC2ContraVoigtStrain(1,1) TLC2ContraVoigtStrain(2,1) TLC2ContraVoigtStrain(3,1)];

% Vector A^{22} (equation (4b))
A22 = [TLC2ContraVoigtStrain(1,2) TLC2ContraVoigtStrain(2,2) TLC2ContraVoigtStrain(3,2)];

% Vector A^{12} (equation (4b))
A12 = .5*[TLC2ContraVoigtStrain(1,3) TLC2ContraVoigtStrain(2,3) TLC2ContraVoigtStrain(3,3)];

%% 4. Compute the derivatives of all entries to the transformation tensor T^{alpha}_{beta},gamma = (A^{alpha}*eLC_{beta}),gamma

% Derivatives T^1_1,gamma (see equations (6a), (6b) in Nitsche - 
% Computation Vol. 7)

% T^1_1,1 (equation (6a))
dT11dxi = dA1Contravariantdxi'*eLC(:,1) + AContravariant(:,1)'*de1LCdxi;

% T^1_1,2 (equation (6b))
dT11deta = dA1Contravariantdeta'*eLC(:,1) + AContravariant(:,1)'*de1LCdeta;

% Derivatives T^1_2,gamma (see equations (7a), (7b) in Nitsche - 
% Computation Vol. 7)

% T^1_2,1 (equation (7a))
dT12dxi = dA1Contravariantdxi'*eLC(:,2) + AContravariant(:,1)'*de2LCdxi;

% T^1_2,1 (equation (7a))
dT12deta = dA1Contravariantdeta'*eLC(:,2) + AContravariant(:,1)'*de2LCdeta;

% Derivatives T^2_1,gamma (see equations (8a), (8b) in Nitsche - 
% Computation Vol. 7)

% T^2_1,1 (equation (8a))
dT21dxi = dA2Contravariantdxi'*eLC(:,1) + AContravariant(:,2)'*de1LCdxi;

% T^2_1,2 (equation (8b))
dT21deta = dA2Contravariantdeta'*eLC(:,1) + AContravariant(:,2)'*de1LCdeta;

% Derivatives T^2_2,gamma (see equations (9a), (9b) in Nitsche - 
% Computation Vol. 7)

% T^2_2,1 (equation (9a))
dT22dxi = dA2Contravariantdxi'*eLC(:,2) + AContravariant(:,2)'*de2LCdxi;

% T^2_2,1 (equation (9a))
dT22deta = dA2Contravariantdeta'*eLC(:,2) + AContravariant(:,2)'*de2LCdeta;

%% 5. Form the matrix dA (see equation (2) in Nitsche - Computation Vol. 7)

% Compute component dA11/dxi (see equation (10) in Nitsche - Computation 
% Vol. 7)
dA11dxi = [2*(AContravariant(:,1)'*eLC(:,1))*dT11dxi 2*(AContravariant(:,1)'*eLC(:,2))*dT12dxi ...
    2*(dT12dxi*(AContravariant(:,1)'*eLC(:,1)) + (AContravariant(:,1)'*eLC(:,2))*dT11dxi)];

% Compute component dA22/deta (see equation (11) in Nitsche - Computation 
% Vol. 7)
dA22deta = [2*(AContravariant(:,2)'*eLC(:,1))*dT21deta 2*(AContravariant(:,2)'*eLC(:,2))*dT22deta ...
    2*(dT21deta*(AContravariant(:,2)'*eLC(:,2)) + (AContravariant(:,2)'*eLC(:,1))*dT22deta)];

% Compute the shear components dA12/dxi and dA12/deta (see equation (12) in
% Nitsche - Computation Vol. 7)

% dA12/dxi
dA12dxi = [dT21dxi*(AContravariant(:,1)'*eLC(:,1)) + (AContravariant(:,2)'*eLC(:,1))*dT11dxi ...
    dT22dxi*(AContravariant(:,1)'*eLC(:,2)) + (AContravariant(:,2)'*eLC(:,2))*dT12dxi ...
    dT21dxi*(AContravariant(:,1)'*eLC(:,2)) + (AContravariant(:,2)'*eLC(:,1))*dT12dxi + ...
    dT22dxi*(AContravariant(:,1)'*eLC(:,1)) + (AContravariant(:,2)'*eLC(:,2))*dT11dxi];

% dA12/deta
dA12deta = [dT21deta*(AContravariant(:,1)'*eLC(:,1)) + (AContravariant(:,2)'*eLC(:,1))*dT11deta ...
    dT22deta*(AContravariant(:,1)'*eLC(:,2)) + (AContravariant(:,2)'*eLC(:,2))*dT12deta ...
    dT21deta*(AContravariant(:,1)'*eLC(:,2)) + (AContravariant(:,2)'*eLC(:,1))*dT12deta + ...
    dT22deta*(AContravariant(:,1)'*eLC(:,1)) + (AContravariant(:,2)'*eLC(:,2))*dT11deta];

% Form matrix dA
dA = [dA11dxi + dA12deta
      dA12dxi + dA22deta];
  
%% 6. Form the derivatives of the the transformation matrix from the contravariant basis to the global Cartesian (see equation (13) in Nitsche - Computation Vol. 7)

% dTLC2ContraVoigtStrain/dxi
dTLC2ContraVoigtStraindxi = ...
    [2*(AContravariant(:,1)'*eLC(:,1))*dT11dxi 2*(AContravariant(:,2)'*eLC(:,1))*dT21dxi 2*(dT11dxi*(AContravariant(:,2)'*eLC(:,1)) + (AContravariant(:,1)'*eLC(:,1))*dT21dxi)
     2*(AContravariant(:,1)'*eLC(:,2))*dT12dxi 2*(AContravariant(:,2)'*eLC(:,2))*dT22dxi 2*(dT12dxi*(AContravariant(:,2)'*eLC(:,2)) + (AContravariant(:,1)'*eLC(:,2))*dT22dxi)
     2*(dT11dxi*(AContravariant(:,1)'*eLC(:,2)) + (AContravariant(:,1)'*eLC(:,1))*dT12dxi) ...
     2*(dT21dxi*(AContravariant(:,2)'*eLC(:,2)) + (AContravariant(:,2)'*eLC(:,1))*dT22dxi) ...
     2*(dT11dxi*(AContravariant(:,2)'*eLC(:,2)) + (AContravariant(:,1)'*eLC(:,1))*dT22dxi + ...
     dT21dxi*(AContravariant(:,1)'*eLC(:,2)) + (AContravariant(:,2)'*eLC(:,1))*dT12dxi)];
 
% dTLC2ContraVoigtStrain/deta
dTLC2ContraVoigtStraindeta = ...
    [2*(AContravariant(:,1)'*eLC(:,1))*dT11deta 2*(AContravariant(:,2)'*eLC(:,1))*dT21deta 2*(dT11deta*(AContravariant(:,2)'*eLC(:,1)) + (AContravariant(:,1)'*eLC(:,1))*dT21deta)
     2*(AContravariant(:,1)'*eLC(:,2))*dT12deta 2*(AContravariant(:,2)'*eLC(:,2))*dT22deta 2*(dT12deta*(AContravariant(:,2)'*eLC(:,2)) + (AContravariant(:,1)'*eLC(:,2))*dT22deta)
     2*(dT11deta*(AContravariant(:,1)'*eLC(:,2)) + (AContravariant(:,1)'*eLC(:,1))*dT12deta) ...
     2*(dT21deta*(AContravariant(:,2)'*eLC(:,2)) + (AContravariant(:,2)'*eLC(:,1))*dT22deta) ...
     2*(dT11deta*(AContravariant(:,2)'*eLC(:,2)) + (AContravariant(:,1)'*eLC(:,1))*dT22deta + ...
     dT21deta*(AContravariant(:,1)'*eLC(:,2)) + (AContravariant(:,2)'*eLC(:,1))*dT12deta)];


%% 7. Compute the Christoffel symbols (see equations (13a)-(13h) in Nitsche - Computation Vol. 9)
Gamma11_1 = dA1Covariant(:,2)'*AContravariant(:,1);
Gamma21_2 = dA1Covariant(:,3)'*AContravariant(:,2);
Gamma22_1 = dA2Covariant(:,2)'*AContravariant(:,1);
Gamma12_1 = dA1Covariant(:,3)'*AContravariant(:,1);
Gamma21_1 = dA1Covariant(:,3)'*AContravariant(:,1);
Gamma22_2 = dA2Covariant(:,2)'*AContravariant(:,2);
Gamma11_2 = dA1Covariant(:,2)'*AContravariant(:,2);
Gamma12_2 = dA1Covariant(:,3)'*AContravariant(:,2);

%% 8. Form the Gamma matrix containing the Christoffel symbols (see equation (12) Nitsche - Computation Vol. 9)
Gamma = [2*Gamma11_1 + Gamma21_2 Gamma22_1               2*Gamma12_1 + Gamma21_1 + Gamma22_2
         Gamma11_2               Gamma12_1 + 2*Gamma22_2 Gamma11_1 + 2*Gamma21_2 + Gamma12_2];

end
