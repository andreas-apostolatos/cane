function [gCovariant, g3Tilde, dA, g3, gabCovariant, dgCovariant, CurvatureVoigt] = ...
    computeMetricsForKirchhoffLoveShellNonLinearCurrent ...
    (i, p, xi, Xi, j, q, eta, Eta, CP, dR, ddR)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns entities of the current geometry needed for nonlinear stiffness 
% matrix of the Kirchhoff-Love shell
%
%          Input :
%            p,q : NURBS polynomial degrees
%            i,j : knot span indices
%         xi,eta : NURBS coodinates of evaluation
%         Xi,Eta : Knot vectors
%             CP : Control Point coordinates of the deformed configuration
%             dR : The first derivatives of the NURBS basis functions for 
%                  the deformed configuration
%            ddR : The second derivatives of the NURBS basis functions for 
%                  the deformed configuration
%
%         Output :
%     gCovariant : The covariant basis of the deformed configuration
%        g3Tilde : The not normalized base vector orthogonal to the surface
%                  for the deformed configuration
%             dA : The differential element area for the deformed
%                  configuration
%             g3 : The normalized normal to the surface base vector for the
%                  deformed configuration
%    dgCovariant : The derivatives of the covariant base vectors
% CurvatureVoigt : The curvature coefficients for the deformed 
%                  configuration in Voigt notation
%
% Function layout :
%
% 1. Compute covariant base vectors g and second derivatives of the position vector
%
% 2. Compute the normal to the surface base vector
%
% 3. Compute the differential element area of the deformed element
%
% 4. Compute the unit normal to the surface base vector
%
% 5. Compute the covariant metric tensor coefficients
%
% 6. Compute the curvature coefficients arranged in Voigt notation
%
%% Function main body

%% 1. Compute covariant base vectors g and second derivatives of the position vector
[gCovariant,dgCovariant] = computeBaseVectorsAndDerivatives2DGivenBasisFunctions(i,p,xi,Xi,j,q,eta,Eta,CP,dR,ddR);

%% 2. Compute the normal to the surface base vector
g3Tilde = cross(gCovariant(:,1),gCovariant(:,2));

%% 3. Compute the differential element area of the deformed element
dA = norm(g3Tilde);

%% 4. Compute the unit normal to the surface base vector
g3 = g3Tilde/dA;

%% 5. Compute the covariant metric tensor coefficients
gabCovariant = zeros(2,2);
gabCovariant(1,1) = gCovariant(:,1)'*gCovariant(:,1);
gabCovariant(1,2) = gCovariant(:,1)'*gCovariant(:,2);
gabCovariant(2,1) = gabCovariant(1,2);
gabCovariant(2,2) = gCovariant(:,2)'*gCovariant(:,2);

%% 6. Compute the curvature coefficients arranged in Voigt notation
CurvatureVoigt = dgCovariant'*g3;

end
