function [BOperatorOmega1, BOperatorOmega2, commonBOperator, Bab] = ...
    computePostprocBOperatorMatrix4RotationsIGAKirchhoffLoveShell ...
    (RMatrix, DRMatrixDxi, DRMatrixDeta, A3Tilde, dA3KommaAlpha, ...
    AContravariant, BVoigt)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the B-operator matrices for the rotations around base vectors A1
% and A2, respectively, at a parametric location of the shell where the 
% basis functions are computed.
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
%         A3Tilde : The not normalized surface normal base vector
%   dA3KommaAlpha : The parametric derivatives of the surface normal vector
%                   = [dA3/dxi dA3/deta]
%  AContravariant : The contravariant base vectors of the curvilinear 
%                   coordinate system
%          BVoigt : The curvature coefficients sorted in vector 
%                   BVoigt = [B11 B22 2*B12]';
%          
%          Output :
% BOperatorOmega1 : The B-operator matrix for the rotation component around
%                   base vector A_1
% BOperatorOmega2 : The B-operator matrix for the rotation component around
%                   base vector A_2
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
% 5. Compute the B-operator for the rotation component around base vector A_1
%
% 6. Compute the B-operator for the rotation component around base vector A_2
%
%% Function main body

%% 1. Compute the transformation matrix of the surface normal vector from the local to the global Cartesian basis
A3 = A3Tilde/norm(A3Tilde);
T3Cov2GC = A3';

%% 2. Compute the derivatives of the transformation matrix of the surface normal vector from the local to the global Cartesian basis w.r.t. both parametric directions, namely, dT3Cov2GC = [d(T3Cov2GC)/d(theta^1) d(T3Cov2GC)/d(theta^2)]
dT3Cov2GC = dA3KommaAlpha';

%% 3. Compute the covariant components of the curvature tensor of the undeformed configuration 
Bab = [BVoigt(1) BVoigt(3)
       BVoigt(3) BVoigt(2)];

%% 4. Compute the B-operator matrix which is common in both the bending and the twisting rotations
commonBOperator = [dT3Cov2GC(1,:)*RMatrix
                   dT3Cov2GC(2,:)*RMatrix] + [T3Cov2GC*DRMatrixDxi
                                              T3Cov2GC*DRMatrixDeta] + Bab'*AContravariant'*RMatrix;
                                        
%% 5. Compute the B-operator for the rotation component around base vector A_1
BOperatorOmega1 = (1/norm(A3Tilde))*[0 1]*commonBOperator;

%% 6. Compute the B-operator for the rotation component around base vector A_2
BOperatorOmega2 = (1/norm(A3Tilde))*[-1 0]*commonBOperator;

end
