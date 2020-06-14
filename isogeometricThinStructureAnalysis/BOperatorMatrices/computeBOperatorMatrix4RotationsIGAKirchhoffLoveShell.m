function [BOperatorOmegaT, BOperatorOmegaN, commonBOperator, Bab] = ...
    computeBOperatorMatrix4RotationsIGAKirchhoffLoveShell ...
    (RMatrix, DRMatrixDxi, DRMatrixDeta, A3, dA3KommaAlpha, AContravariant, ...
    BVoigt, nCovariant, tCovariant)
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
% are computed and along the curve on the shell which is uniquely described
% by the tangent and normal vectors TCovariant, NCovariant respectively.
%
%           Input :
%         RMatrix : Matrix containing the basis functions corresponding to
%                   a 3D displacement field
%     DRMatrixDxi : Matrix containing the derivatives of the basis 
%                   functions corresponding to a 3D displacement field with
%                   respect to parameter xi
%    DRMatrixDeta : Matrix containing the derivatives of the basis 
%                   functions corresponding to a 3D displacement field with
%                   respect to parameter eta
%              A3 : The unit normal to the NURBS surface base vector
%   dA3KommaAlpha : The parametric derivatives of the surface normal vector
%                   = [dA3/dxi dA3/deta]
%  AContravariant : The contravariant base vectors of the curvilinear 
%                   coordinate system
%          BVoigt : The curvature coefficients sorted in vector 
%                   BVoigt = [B11 B22 2*B12]';
%      nCovariant : The normal to the shell boundary vector in the
%                   covariant basis
%      tCovariant : The tangent to the shell boundary vector in the
%                   covariant basis
%          
%          Output :
% BOperatorOmegaT : The B-operator matrix for the bending rotations
% BOperatorOmegaN : The B-operator matrix for the twisting rotations
% commonBOperator : The common B-operator matrix for both the bending and
%                   the twisting rotations
%             Bab : The curvature tensor in the contravariant basis
%                   arranged in a square matrix
%
% Function layout :
%
% 1. Compute the transformation matrix of the surface normal vector from the local to the global Cartesian basis
%
% 2. Compute the derivatives of the transformation matrix of the surface normal vector from the local to the global Cartesian basis w.r.t. both parametric directions, namely, dT3Cov2GC = [d(T3Cov2GC)/d(theta^1) d(T3Cov2GC)/d(theta^2)]
% 
% 3. Compute the covariant components of the curvature tensor of the undeformed configuration
%
% 4. Compute the B-operator matrix which is common in both the bending and the twisting rotations
%
% 5. Compute the B-operator for the bending rotations
%
% 6. Compute the B-operator for the twisting rotations
%
%% Function main body

%% 1. Compute the transformation matrix of the surface normal vector from the local to the global Cartesian basis
T3Cov2GC = A3';

%% 2. Compute the derivatives of the transformation matrix of the surface normal vector from the local to the global Cartesian basis w.r.t. both parametric directions, namely, dT3Cov2GC = [d(T3Cov2GC)/d(theta^1) d(T3Cov2GC)/d(theta^2)]
dT3Cov2GC = dA3KommaAlpha';

%% 3. Compute the covariant components of the curvature tensor of the undeformed configuration 
Bab = [BVoigt(1) BVoigt(3)
       BVoigt(3) BVoigt(2)];

%% 4. Compute the B-operator matrix which is common in both the bending and the twisting rotations
%
% commonBOperator1 = dT3Cov2GC*RMatrix;
%
% commonBOperator2 = [T3Cov2GC*DRMatrixDxi
%                     T3Cov2GC*DRMatrixDeta];
%
% commonBOperator3 = Bab'*AContravariant'*RMatrix;
%
% commonBOperator = commonBOperator1 + commonBOperator2 + commonBOperator3;
%
commonBOperator = [dT3Cov2GC(1,:)*RMatrix
                   dT3Cov2GC(2,:)*RMatrix] + [T3Cov2GC*DRMatrixDxi
                                              T3Cov2GC*DRMatrixDeta] + Bab'*AContravariant'*RMatrix;
                                          
%% 5. Compute the B-operator for the bending rotations
BOperatorOmegaT = - nCovariant'*commonBOperator;

%% 6. Compute the B-operator for the twisting rotations
BOperatorOmegaN = tCovariant'*commonBOperator;

end
