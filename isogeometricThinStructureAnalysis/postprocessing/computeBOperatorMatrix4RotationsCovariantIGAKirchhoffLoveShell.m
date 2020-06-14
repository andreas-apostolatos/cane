function BOperatorOmegaLC = ...
    computeBOperatorMatrix4RotationsCovariantIGAKirchhoffLoveShell ...
    (p, q, dR, A3, dA3, dA, ACovariant, dACovariant)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the B-operator matrices for the bending and the twisting
% rotations at a parametric location of the shell where the basis functions
% are computed and along the curve on the shell in both parametric
% directions over the local Cartesian basis.
%
%            Input :
%              p,q : The polynomial degrees of the B-Spline surface in both
%                    parametric directions
%               dR : The basis functions and up to their first derivatives    
%               A3 : The unit normal to the NURBS surface base vector
%              dA3 : The derivatives of the surface normal vector
%                    = [dA3/dxi dA3/deta]
%               dA : The differential element area
%       ACovariant : the covariant base vectors
%      dACovariant : The parametric derivatives of the covariant base
%                    vectors = [dA1/dxi dA2/deta dA1/deta(=dA2/dxi)] 
%          
%           Output :
% BOperatorOmegaLC : The B-operator matrix for the rotation vector in the
%                    local Cartesian basis
%
% Function layout :
%
% 0. Read input
%
% 1. Compute the matrices containing the basis functions and their derivatives
%
% 2. Compute the covariant metric coefficients
%
% 3. Compute the contravariant base vectors
%
% 4. Compute the Voigt curvature tensor in the contravariant basis [B11 B22 2*B12]'
%
% 5. Compute the Voigt curvature tensor in the contravariant basis [B11 B22 2*B12]'
%
% 6. Compute the transformation matrix of the surface normal vector from the local to the global Cartesian basis
%
% 7. Compute the derivatives of the transformation matrix of the surface normal vector from the local to the global Cartesian basis w.r.t. both parametric directions, namely, dT3Cov2GC = [d(T3Cov2GC)/d(theta^1) d(T3Cov2GC)/d(theta^2)]
%
% 8. Compute the covariant components of the curvature tensor of the undeformed configuration 
%
% 9. Compute the B-operator matrix which is common in the rotations in both parametric directions
%
% 10. Compute the B-operator for the roations in the xi-direction
%
% 11. Compute the B-operator for the roations in the eta-direction
%
% 12. Form the complete B-operator matrix for the rotations in the covariant basis
%
% 13. Form the complete B-operator matrix for the rotations in the local Cartesian basis
%
%% Function main body

%% 0. Read input

% Number of local Control Points
noCPsLoc = (p + 1)*(q + 1);

% Number of local DOFs
noDOFsLoc = 3*noCPsLoc;

% Initialize auxiliary arrays
RMatrix = zeros(3,noDOFsLoc);
DRMatrixDxi = zeros(3,noDOFsLoc);
DRMatrixDeta = zeros(3,noDOFsLoc);

%% 1. Compute the matrices containing the basis functions and their derivatives
for counterBasis = 1:noCPsLoc
    % Matrix containing the basis functions
    RMatrix(1,3*counterBasis-2) = dR(counterBasis,1);
    RMatrix(2,3*counterBasis-1) = dR(counterBasis,1);
    RMatrix(3,3*counterBasis) = dR(counterBasis,1);

    % Matrix containing the derivatives of the basis functions
    % wrt xi-direction
    DRMatrixDxi(1,3*counterBasis-2) = dR(counterBasis,2);
    DRMatrixDxi(2,3*counterBasis-1) = dR(counterBasis,2);
    DRMatrixDxi(3,3*counterBasis) = dR(counterBasis,2);

    % Matrix containing the derivatives of the basis functions
    % wrt eta-direction
    DRMatrixDeta(1,3*counterBasis-2) = dR(counterBasis,4);
    DRMatrixDeta(2,3*counterBasis-1) = dR(counterBasis,4);
    DRMatrixDeta(3,3*counterBasis) = dR(counterBasis,4);
end

%% 2. Compute the covariant metric coefficients
AabCov = [ACovariant(:,1) ACovariant(:,2)]'*[ACovariant(:,1) ACovariant(:,2)];

%% 3. Compute the contravariant base vectors
AContravariant = (AabCov\[ACovariant(:,1) ACovariant(:,2)]')';

%% 4. Compute the local Cartesian basis
eLC = computeLocalCartesianBasis4BSplineSurface(ACovariant,AContravariant);

%% 5. Compute the Voigt curvature tensor in the contravariant basis [B11 B22 2*B12]'
BVoigt = [dACovariant(:,1) dACovariant(:,2) dACovariant(:,3)]'*A3;

%% 6. Compute the transformation matrix of the surface normal vector from the local to the global Cartesian basis
T3Cov2GC = A3';

%% 7. Compute the derivatives of the transformation matrix of the surface normal vector from the local to the global Cartesian basis w.r.t. both parametric directions, namely, dT3Cov2GC = [d(T3Cov2GC)/d(theta^1) d(T3Cov2GC)/d(theta^2)]
dT3Cov2GC = dA3';

%% 8. Compute the covariant components of the curvature tensor of the undeformed configuration 
Bab = [BVoigt(1) BVoigt(3)
       BVoigt(3) BVoigt(2)];

%% 9. Compute the B-operator matrix which is common in the rotations in both parametric directions
commonBOperator = [dT3Cov2GC(1,:)*RMatrix
                   dT3Cov2GC(2,:)*RMatrix] + [T3Cov2GC*DRMatrixDxi
                                              T3Cov2GC*DRMatrixDeta] + Bab'*AContravariant'*RMatrix;
                                        
%% 10. Compute the B-operator for the roations in the xi-direction
permutationVct1 = - 1/dA*[0
                          1];
BOperatorOmega1 = - permutationVct1'*commonBOperator;

%% 11. Compute the B-operator for the roations in the eta-direction
permutationVct2 = 1/dA*[1
                        0];
BOperatorOmega2 = - permutationVct2'*commonBOperator;

%% 12. Form the complete B-operator matrix for the rotations in the covariant basis
BOperatorOmegaCov = [BOperatorOmega1
                     BOperatorOmega2];
              
%% 13. Form the complete B-operator matrix for the rotations in the local Cartesian basis
BOperatorOmegaLC = [ACovariant(:,1)'*eLC(:,1) ACovariant(:,2)'*eLC(:,1)
                    ACovariant(:,1)'*eLC(:,2) ACovariant(:,2)'*eLC(:,2)]*BOperatorOmegaCov;

end
