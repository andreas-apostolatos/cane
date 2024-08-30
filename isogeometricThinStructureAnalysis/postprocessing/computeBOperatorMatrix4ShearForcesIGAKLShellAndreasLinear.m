function BShearForces = ...
    computeBOperatorMatrix4ShearForcesIGAKLShellAndreasLinear...
    (p, q, dR, dACovariant1, dACovariant2, A3Tilde, Db)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the B-operator for the shear forces in the local Cartesian basis 
% corresponding to the linear formulation of the Kirchhoff-Love shell in 
% isogeometric analysis derived by Andreas.
%
%           Input : 
%             p,q : polynomial degrees
%              dR : The basis function and up to its third mixed
%                   derivatives
%	 dACovariant1 : The base vector A1 and up to its second mixed 
%                   derivatives :
%                   = [A1 dA1/dxi d^2A1/dxi^2 dA1/deta(=dA2/dxi) 
%                      d^2A1/detadxi(=d^2A2/dxi^2) 
%                      d^2A1/deta^2(=d^2A2/dxideta)]
%    dACovariant2 : The base vector A1 and up to its second mixed 
%                   derivatives (symmetry with respect to the A1 base 
%                   vector is taken into consideration)
%                   = [A2 dA2/deta d^2A2/deta^2]
%         A3Tilde : The not normalized surface normal
%              Db : The material matrix of the bending part of the
%                   stiffness
%
%          Output :
%    BShearForces : The B-operator matrix for the shear forces 2 x nDOFsEl
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the basis functions
% ->
%    1i. Form matrix dRdxi
%
%   1ii. Form matrix dRdeta
%
%  1iii. Form matrix ddRdxidxi
%
%   1iv. Form matrix ddRdetadeta
%
%    1v. Form matrix ddRdxideta
%
%   1vi. Form matrix dddRdxidxidxi
%
%  1vii. Form matrix dddRdxidxideta
%
% 1viii. Form matrix dddRdetadetadxi
%
%   1ix. Form matrix dddRdetadetadeta
% <-
%
% 2. Compute the covariant metric coefficients
%
% 3. Compute the contravariant base vectors
%
% 4. Compute the local Cartesian basis
%
% 5. Compute the transformation matrices from the local Cartesian to the contravariant bases
%
% 6. Compute the derivatives of the surface normal vector and the element surface area
%
% 7. Compute the curvature coefficients B_{alpha beta} = A_{alpha,beta}*A_{3}
%
% 8. Compute the B-operator matrices for the strain tensors, the change in the curvature tensors and the derivatives of the bending strain tensors in the contravariant basis
%
% 9. Compute the derivatives of the tranformation matrices needed for the computation of the B-operator to the shear force tensor
%
% 10. Compute the B-operator matrices for the shear forces [q^1 q^2] in the covariant basis
%
% 11. Compute the transformation matrix from the covariant to the local Cartesian basis for the shear forces
%
% 12. Transform the B-operator matrix to the local Cartesian basis
% 
%% Function main body

%% 0. Read input

% Compute the surface normal vector
A3 = A3Tilde/norm(A3Tilde);

% The number of Control Points
nCPsEl = (p + 1)*(q + 1);

% The number of Control Points
nDOFsEl = 3*nCPsEl;

% Initialize auxiliary arrays

% RKomma1
dRdxi = zeros(3,nDOFsEl);

% RKomma2
dRdeta = zeros(3,nDOFsEl);

% RKomma1Komma1
ddRdxidxi = zeros(3,nDOFsEl);

% RKomma2Komma2
ddRdetadeta = zeros(3,nDOFsEl);

% RKomma1Komma2s
ddRdxideta = zeros(3,nDOFsEl);

% RKomma1Komma1Komma1
dddRdxidxidxi = zeros(3,nDOFsEl);

% RKomma1Komma1Komma2
dddRdxidxideta = zeros(3,nDOFsEl);

% RKomma2Komma2Komma1
dddRdetadetadxi = zeros(3,nDOFsEl);

% RKomma2Komma2Komma2
dddRdetadetadeta = zeros(3,nDOFsEl);

%% 1. Loop over all the basis functions
for i = 1:nCPsEl
    %% 1i. Form matrix dRdxi
    dRdxi(1,3*i-2) = dR(i,2);
    dRdxi(2,3*i-1) = dR(i,2);
    dRdxi(3,3*i) = dR(i,2);

    %% 1ii. Form matrix dRdeta
    dRdeta(1,3*i-2) = dR(i,5);
    dRdeta(2,3*i-1) = dR(i,5);
    dRdeta(3,3*i) = dR(i,5);

    %% 1iii. Form matrix ddRdxidxi
    ddRdxidxi(1,3*i-2) = dR(i,3);
    ddRdxidxi(2,3*i-1) = dR(i,3);
    ddRdxidxi(3,3*i) = dR(i,3);

    %% 1iv. Form matrix ddRdetadeta
    ddRdetadeta(1,3*i-2) = dR(i,8);
    ddRdetadeta(2,3*i-1) = dR(i,8);
    ddRdetadeta(3,3*i) = dR(i,8);

    %% 1v. Form matrix ddRdxideta
    ddRdxideta(1,3*i-2) = dR(i,6);
    ddRdxideta(2,3*i-1) = dR(i,6);
    ddRdxideta(3,3*i) = dR(i,6);

    %% 1vi. Form matrix dddRdxidxidxi
    dddRdxidxidxi(1,3*i-2) = dR(i,4);
    dddRdxidxidxi(2,3*i-1) = dR(i,4);
    dddRdxidxidxi(3,3*i) = dR(i,4);

    %% 1vii. Form matrix dddRdxidxideta
    dddRdxidxideta(1,3*i-2) = dR(i,7);
    dddRdxidxideta(2,3*i-1) = dR(i,7);
    dddRdxidxideta(3,3*i) = dR(i,7);

    %% 1viii. Form matrix dddRdetadetadxi
    dddRdetadetadxi(1,3*i-2) = dR(i,9);
    dddRdetadetadxi(2,3*i-1) = dR(i,9);
    dddRdetadetadxi(3,3*i) = dR(i,9);

    %% 1ix. Form matrix dddRdetadetadeta
    dddRdetadetadeta(1,3*i-2) = dR(i,10);
    dddRdetadetadeta(2,3*i-1) = dR(i,10);
    dddRdetadetadeta(3,3*i) = dR(i,10);
end

%% 2. Compute the covariant metric coefficients
AabCov = [dACovariant1(:,1) dACovariant2(:,1)]'*[dACovariant1(:,1) dACovariant2(:,1)];

%% 3. Compute the contravariant base vectors
AContravariant = (AabCov\[dACovariant1(:,1) dACovariant2(:,1)]')';

%% 4. Compute the local Cartesian basis
eLC = computeLocalCartesianBasis4BSplineSurface([dACovariant1(:,1) dACovariant2(:,1)],AContravariant);

%% 5. Compute the transformation matrices from the local Cartesian to the contravariant bases
TLC2ContraVoigt = computeTFromContra2LocalCartesian4VoigtStrainIGAKLShell(eLC,AContravariant);

%% 6. Compute the derivatives of the surface normal vector and the element surface area
[dA3,dAKommaAlpha] = computeParametricDrvsSurfaceNormalOnBSplineSurface...
    ([dACovariant1(:,1) dACovariant2(:,1)],[dACovariant1(:,2) dACovariant2(:,2) dACovariant1(:,4)],...
    A3,norm(A3Tilde));

%% 7. Compute the curvature coefficients B_{alpha beta} = A_{alpha,beta}*A_{3}
BV = [dACovariant1(:,2) dACovariant2(:,2) dACovariant1(:,4)]'*A3;

%% 8. Compute the B-operator matrices for the strain tensors, the change in the curvature tensors and the derivatives of the bending strain tensors in the contravariant basis
[~,BBendingStrainContra,dBBendingStrainContradxi,dBBendingStrainContradeta] = ...
                computeBOperatorMatrices4IGAKirchhoffLoveShellLinear...
                (dRdxi,dRdeta,ddRdxidxi,ddRdetadeta,ddRdxideta,dddRdxidxidxi,...
                 dddRdxidxideta,dddRdetadetadxi,dddRdetadetadeta,dACovariant1,...
                 dACovariant2,A3Tilde,dA3(:,1),dA3(:,2),dAKommaAlpha,BV);
             
%% 9. Compute the derivatives of the tranformation matrices needed for the computation of the B-operator to the shear force tensor
[A11,A12,A22,dA,dTLC2ContraVoigtStraindxi,dTLC2ContraVoigtStraindeta] = ...
                computeDerivativesTransformationMtcesDDMNitscheIGAKLShell...
                ([dACovariant1(:,1) dACovariant1(:,2) dACovariant1(:,4)],[dACovariant2(:,1) dACovariant2(:,2)],...
                eLC,AContravariant,TLC2ContraVoigt);
            
%% 10. Compute the B-operator matrices for the shear forces [q^1 q^2] in the covariant basis
BShearForcesCovariant = (dA*Db*TLC2ContraVoigt + ...
                            [A11
                             A12]*Db*dTLC2ContraVoigtStraindxi + ...
                            [A12
                             A22]*Db*dTLC2ContraVoigtStraindeta)*BBendingStrainContra + ...
                             [A11
                              A12]*Db*TLC2ContraVoigt*dBBendingStrainContradxi + ...
                              [A12
                               A22]*Db*TLC2ContraVoigt*dBBendingStrainContradeta;
               
%% 11. Compute the transformation matrix from the covariant to the local Cartesian basis for the shear forces
TCov2LC = [dACovariant1(:,1)'*eLC(:,1) dACovariant2(:,1)'*eLC(:,1)
           dACovariant1(:,1)'*eLC(:,2) dACovariant2(:,1)'*eLC(:,2)];

%% 12. Transform the B-operator matrix to the local Cartesian basis
BShearForces = TCov2LC*BShearForcesCovariant;

end
