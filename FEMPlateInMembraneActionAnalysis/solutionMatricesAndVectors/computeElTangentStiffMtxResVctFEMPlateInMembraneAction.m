function [tanMtxMatEl, tanMtxGeoEl, resIntEl, FBodyEl] = ...
    computeElTangentStiffMtxResVctFEMPlateInMembraneAction ...
    (uEl, dN, bF, propParameters, C, DetJxxi, GW)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the material, the geometric, the mass matrix, the internal 
% residual and external body force vector at the Gauss point corresponding
% to the geometrically nonlinear plane stress formulation.
%
%          Input :
%            uEl : The displacement vector at the element level arranged as
%                  uEl = [u1x u1y; u2x u2y; u3x u3y]
%             dN : The basis functions and their derivatives
%             bF : The body force vector
% propParameters : The parameters of the problem
%                    .E : The Young's modulus
%                  .nue : The Poisson's ratio
%                  .rho : The material density
%              C : The material matrix of the plate in membrane action 
%                  problem
%        DetJxxi : The element area at the Gauss Point
%             GW : The corresponding Gauss weight
% 
%         Output :
%    tanMtxMatEl : The material stiffness at the Gauss point
%    tanMtxGeoEl : The geometric stiffness at the Gauss point
%       resIntEl : The internal residual K(u)*u at the Gauss point
%        FBodyEl : The body force vector at the Gauss point
%
% Function layout :
%
% 1. Compute the deformation gradient tensor F at the Gauss point
%
% 2. Compute the B-operator matrix corresponding to the current configuration at the Gauss point
%
% 3. Compute the local material stiffness matrix at the Gauss point
%
% 4. Compute the Green-Lagrange strain tensor and re-arrange it to the Voigt notation at the Gauss point
%
% 5. Compute the Cauchy stress tensor out of the Voigt Green-Lagrange strain vector at the Gauss point
%
% 6. Compute the element geometric stiffness matrix at the Gauss point
%
% 7. Compute the internal residual K(u)*u at the Gauss point
%
% 8. Compute the body force vector at the Gauss point
%
%% Function main body

%% 1. Compute the deformation gradient tensor F at the Gauss point
        
% Re-arrange the derivatives of the basis functions as
%
%      |dN1/dX dN1/dy|
% dN = |dN2/dX dN2/dy|
%      |dN3/dX dN3/dy|
%
dNdX = [dN(1,2:3)
        dN(2,2:3)
        dN(3,2:3)];

% Compute the deformation gradient tensor
%
%      |1 0|   |dx/dX    dx/dY|
%  F = |0 1| + |dy/dX    dy/dY| 
%
FDefGrad = eye(2) + uEl'*dNdX;

%% 2. Compute the B-operator matrix corresponding to the current configuration at the Gauss point
B = [dNdX(1,1)*FDefGrad(1,1)                           dNdX(1,1)*FDefGrad(2,1)                           dNdX(2,1)*FDefGrad(1,1)                           dNdX(2,1)*FDefGrad(2,1)                           dNdX(3,1)*FDefGrad(1,1)                           dNdX(3,1)*FDefGrad(2,1)
     dNdX(1,2)*FDefGrad(1,2)                           dNdX(1,2)*FDefGrad(2,2)                           dNdX(2,2)*FDefGrad(1,2)                           dNdX(2,2)*FDefGrad(2,2)                           dNdX(3,2)*FDefGrad(1,2)                           dNdX(3,2)*FDefGrad(2,2)
     dNdX(1,1)*FDefGrad(1,2) + dNdX(1,2)*FDefGrad(1,1) dNdX(1,1)*FDefGrad(2,2) + dNdX(1,2)*FDefGrad(2,1) dNdX(2,1)*FDefGrad(1,2) + dNdX(2,2)*FDefGrad(1,1) dNdX(2,1)*FDefGrad(2,2) + dNdX(2,2)*FDefGrad(2,1) dNdX(3,1)*FDefGrad(1,2) + dNdX(3,2)*FDefGrad(1,1) dNdX(3,1)*FDefGrad(2,2) + dNdX(3,2)*FDefGrad(2,1)];

%% 3. Compute the local material stiffness matrix at the Gauss point
tanMtxMatEl = B'*C*B*DetJxxi*GW;

%% 4. Compute the Green-Lagrange strain tensor and re-arrange it to the Voigt notation at the Gauss point
epsilonGLTensor = 0.5*(FDefGrad'*FDefGrad - eye(2));
epsilonGLVoigt = [epsilonGLTensor(1,1)
                  epsilonGLTensor(2,2)
                  2*epsilonGLTensor(1,2)];

%% 5. Compute the Cauchy stress tensor out of the Voigt Green-Lagrange strain vector at the Gauss point
stressCauchyVoigt = C*epsilonGLVoigt;
stressCauchyTensor = [stressCauchyVoigt(1)  stressCauchyVoigt(3)
                      stressCauchyVoigt(3)  stressCauchyVoigt(2)];

%% 6. Compute the element geometric stiffness matrix at the Gauss point
H = dNdX*stressCauchyTensor*dNdX';
tanMtxGeoEl = [H(1,1) 0      H(1,2)  0      H(1,3)  0
               0      H(1,1) 0       H(1,2) 0       H(1,3)
               H(2,1) 0      H(2,2)  0      H(2,3)  0
               0      H(2,1) 0       H(2,2) 0       H(2,3)
               H(3,1) 0      H(3,2)  0      H(3,3)  0
               0      H(3,1) 0       H(3,2) 0       H(3,3)]*DetJxxi*GW;

%% 7. Compute the internal residual K(u)*u at the Gauss point
resIntEl = B'*stressCauchyVoigt*DetJxxi*GW;

%% 8. Compute the body force vector at the Gauss point
NMtx = [dN(1,1) 0       dN(2,1) 0       dN(3,1) 0
        0       dN(1,1) 0       dN(2,1) 0       dN(3,1)];
% FBodyEl = NMtx'*bF(1:2,1)*DetJxxi*GW;
FBodyEl = NMtx'*squeeze(bF(1,1,1:2))*DetJxxi*GW;

end
