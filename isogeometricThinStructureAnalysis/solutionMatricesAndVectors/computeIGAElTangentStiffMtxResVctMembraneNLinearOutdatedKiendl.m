function [tanMtxEl, resVctEl] = ...
    computeIGAElTangentStiffMtxResVctMembraneNLinearOutdatedKiendl ...
    (i, p, xi, Xi, j, q, eta, Eta, CP, dR, ACovariant, aCovariant, ...
    thickness, Dm, prestress)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the symmetric part of the tangent stiffness matrix (only the 
% upper triangular part) and load vector corresponding to the isogeometric 
% membrane element. The implementation style is as in the PhD thesis of 
% Josef Kield and inefficient but very clear to read.
%
%                 Input :
%                   i,j : The knot span indices
%                xi,eta : The parametric coordinates where the basis
%                         functions are computed
%                Xi,Eta : The knot vector in xi- and eta-directions
%                    CP : The array of Control Point coordinates and
%                         weights
%                   p,q : The polynomial degrees of the B-Spline patch
%                    dR : The B-Spline basis functions and their first 
%                         derivatives = [R dR/dxi dR/deta]
%            ACovariant : The covariant base vectors of the reference 
%                         configuration
%            aCovariant : The covariant base vectors of the current 
%                         configuration
%             thickness : The thickness of the membrane
%                    Dm : The membrane material stiffness matrix
%             prestress : On the actual prestress of the membrane
%                          .computeParametricCoordinates : Function handle
%                                                          to the 
%                                                          computation of 
%                                                          the coordinates
%                                                          of the system
%                                                          where the
%                                                          prestress values
%                                                          are prescribed
%                                    .computeBaseVectors : Function handle 
%                                                          to the 
%                                                          computation of 
%                                                          the base vectors 
%                                                          of the system 
%                                                          where the 
%                                                          prestress values 
%                                                          are prescribed
%                                           .voigtVector : Vector 
%                                                          containing the 
%                                                          prestress 
%                                                          coefficients in
%                                                          Voigt vector on
%                                                          the system
%                                                          defined by the
%                                                          above function
%                                                          handles or if
%                                                          those handles
%                                                          are not defined
%                                                          over the local
%                                                          Cartesian system
%
%                Output :
%              tanMtxEl : Element tangent stiffness matrix
%              resVctEl : Element internal residual load vector
%
% Function layout :
%
% 0. Read input
%
% 1. Compute metrics for the reference and the current configuration
%
% 2. Compute the prestress values on the local Cartesian coordinate system
%
% 3. Compute the nonlinear Green-Lagrange membrane strain in the contravariant basis
%
% 4. Transform the nonlinear Green-Lagrange membrane strain in the local Cartesian basis
%
% 5. Loop over the DOFs to compute the first variation wrt DOFs
% ->
%    5i. Get local node number kr and dof direction dirr
%
%   5ii. Compute first variation of the Green-Lagrange strain wrt the DOFS in the contravariant system
%
%  5iii. Compute first variation of the Green-Lagrange strain wrt the DOFS in the contravariant system
% <-
%
% 6. Loop over the r-DOFs to compute the second variation wrt DOFs
% ->
%    6i. Get local node number kr and dof direction dirr for DoF ur
%
%   6ii. Loop over the s-DOFs to compute the second variation wrt DOFs
%   ->
%        6ii.1. Get local node number ks and dof direction dirs for DoF us
%
%        6ii.2. Compute second variation of the Green-Lagrange strain wrt DOFs in the contravariant system
%
%        6ii.3. Transform the second variation of the Green-Lagrange strain wrt DOFs onto the local Cartesian basis
%   <-
% <-
%
% 7. Loop over the r-DOFs to compute the second variation of the 2nd Piola-Kirchhoff stress wrt the DOFs
% ->
%    7i. Compute the 2nd Piola-Kirchhoff stress in the local Cartesian basis
%
%   7ii. Compute the first variation of the 2nd Piola-Kirchhoff stress wrt the DOFs in the local Cartesian basis
%
%  7iii. Compute the tangent stiffness matrix
%
%   7iv. Compute the residual vector
% <-
% 
%% 0. Read input

% Number of DoFs at the element level
numDOFsEl = 3*(p + 1)*(q + 1);

% Initialize element tangent stiffness matrix
tanMtxEl = zeros(numDOFsEl, numDOFsEl);

% Initialize element internal residual load vector
resVctEl = zeros(numDOFsEl, 1);

% Initialize arrays on the first variation of the strain and the curvature
% w.r.t. the DoFs
dEpsilonLC = zeros(3, numDOFsEl);

% Initialize arrays on the second variation of the strain wrt the DoFs
ddEpsilonLC = zeros(3, numDOFsEl, numDOFsEl);

%% 1. Compute metrics for the reference and the current configuration

% For the reference configuration :
% _________________________________

% Compute the covariant metric coefficients
AabCovariant = zeros(2, 2);
AabCovariant(1, 1) = ACovariant(:, 1)'*ACovariant(:, 1);
AabCovariant(1, 2) = ACovariant(:, 1)'*ACovariant(:, 2);
AabCovariant(2, 1) = AabCovariant(1, 2);
AabCovariant(2, 2) = ACovariant(:, 2)'*ACovariant(:, 2);

% Compute the contravariant basis
AContravariant = AabCovariant\ACovariant';
AContravariant = AContravariant';

% Compute the local Cartesian basis
eLC = computeLocalCartesianBasis4BSplineSurface ...
    (ACovariant, AContravariant);

% Compute the transformation matrix from the contravariant basis to the
% local Cartesian one
TFromContraToLC = computeTFromContra2LocalCartesian4VoigtStrainIGAKLShell ...
    (eLC, AContravariant);

% For the current configuration :
% _______________________________

% Compute the covariant metric coefficients
aabCovariant = zeros(2, 2);
aabCovariant(1, 1) = aCovariant(:, 1)'*aCovariant(:, 1);
aabCovariant(1, 2) = aCovariant(:, 1)'*aCovariant(:, 2);
aabCovariant(2, 1) = aabCovariant(1, 2);
aabCovariant(2, 2) = aCovariant(:, 2)'*aCovariant(:, 2);

%% 2. Compute the prestress values on the local Cartesian coordinate system
isPrestressOverDefinedSystem = false;
if (isfield(prestress, 'computeParametricCoordinates') && ~isfield(prestress, 'computeBaseVectors')) || ...
        (~isfield(prestress, 'computeParametricCoordinates') && isfield(prestress, 'computeBaseVectors'))
    error('Both or none of function handles parameters.prestress.computeParametricCoordinates and parameters.prestress.computeBaseVectors have to be defined');
elseif isfield(prestress, 'computeParametricCoordinates') && isfield(prestress, 'computeBaseVectors')
    isPrestressOverDefinedSystem = true;
end
if isPrestressOverDefinedSystem
    X = computeCartesianCoordinatesOfAPointOnBSplineSurface ...
        (i, p, xi, Xi, j, q, eta, Eta, CP, dR(:,1));
    theta = prestress.computeParametricCoordinates(X);
    prestressBaseVct = prestress.computeBaseVectors(theta(1, 1), theta(2, 1));
    T2LC = computeT2LocalCartesianBasis(prestressBaseVct, eLC);
    pTilde = T2LC*prestress.voigtVector;
else
    pTilde = prestress.voigtVector;
end

%% 3. Compute the nonlinear Green-Lagrange membrane strain in the contravariant basis
EpsilonContra(1, 1) = 1/2*(aabCovariant(1, 1) - AabCovariant(1, 1));
EpsilonContra(2, 1) = 1/2*(aabCovariant(2, 2) - AabCovariant(2, 2));
EpsilonContra(3, 1) = 1/2*(aabCovariant(1, 2) - AabCovariant(1, 2));

%% 4. Transform the nonlinear Green-Lagrange membrane strain in the local Cartesian basis
EpsilonLC = TFromContraToLC*EpsilonContra;

%% 5. Loop over the DOFs to compute the first variation wrt DOFs
for rDOFs = 1:numDOFsEl
    %% 5i. Get local node number kr and dof direction dirr
    kr = ceil(rDOFs/3);
    dirr = rDOFs - 3*(kr - 1);
  
    %% 5ii. Compute first variation of the Green-Lagrange strain wrt the DOFS in the contravariant system
    dEpsilonContra(1, 1) = dR(kr, 2)*aCovariant(dirr, 1);
    dEpsilonContra(2, 1) = dR(kr, 3)*aCovariant(dirr, 2);
    dEpsilonContra(3, 1) = 1/2*(dR(kr, 2)*aCovariant(dirr, 2) + ...
        aCovariant(dirr, 1)*dR(kr, 3));
    
    %% 5iii. Compute first variation of the Green-Lagrange strain wrt the DOFS in the contravariant system
    dEpsilonLC(:, rDOFs) = TFromContraToLC*dEpsilonContra;
end

%% 6. Loop over the r-DOFs to compute the second variation wrt DOFs
for rDOFs = 1:numDOFsEl
    %% 6i. Get local node number kr and dof direction dirr for DoF ur
    kr = ceil(rDOFs/3);
    dirr = rDOFs - 3*(kr - 1);
    
    %% 6ii. Loop over the s-DOFs to compute the second variation wrt DOFs
    for sDOFs = 1:rDOFs
        %% 6ii.1. Get local node number ks and dof direction dirs for DoF us
        ks = ceil(sDOFs/3);
        dirs = sDOFs - 3*(ks - 1);
    
        %% 6ii.2. Compute second variation of the Green-Lagrange strain wrt DOFs in the contravariant system
        ddEpsilonContra = zeros(3, 1);
        if dirr == dirs
            ddEpsilonContra(1, 1) = dR(kr, 2)*dR(ks, 2);
            ddEpsilonContra(2, 1) = dR(kr, 3)*dR(ks, 3);
            ddEpsilonContra(3, 1) = 1/2*(dR(kr, 2)*dR(ks, 3) + dR(kr, 3)*dR(ks, 2));
        end

        %% 6ii.3. Transform the second variation of the Green-Lagrange strain wrt DOFs onto the local Cartesian basis
        ddEpsilonLC(:, rDOFs, sDOFs) = TFromContraToLC*ddEpsilonContra;
    end
end

%% 7. Loop over the r-DOFs to compute the second variation of the 2nd Piola-Kirchhoff stress wrt the DOFs
for rDOFs = 1:numDOFsEl
    %% 7i. Compute the 2nd Piola-Kirchhoff stress in the local Cartesian basis
    NLC = thickness*pTilde +  Dm*EpsilonLC;

    %% 7ii. Compute the first variation of the 2nd Piola-Kirchhoff stress wrt the DOFs in the local Cartesian basis
    dNLC = Dm*dEpsilonLC(:, rDOFs);
    
    %% 7iii. Compute the tangent stiffness matrix
    for sDOFs = 1:rDOFs
        tanMtxEl(rDOFs, sDOFs) = dNLC'*dEpsilonLC(:, sDOFs) + NLC'*ddEpsilonLC(:, rDOFs, sDOFs);
    end
    
    %% 7iv. Compute the residual vector
    resVctEl(rDOFs, 1) = NLC'*dEpsilonLC(:, rDOFs);
end

end
