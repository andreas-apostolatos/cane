function Ke = computeElStiffMtxHSDTShellLinear ...
    (p, q, dR, GCovariant, dGCovariant, G3Tilde, Dm, Db, Ds, shearCorrection)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    [Your Name] (based on Andreas Apostolatos)
%
%% Function documentation
%
% Returns the element stiffness matrix corresponding to the linear
% isogeometric Higher-order Shear Deformation Theory (HSDT) shell problem
% with 5 DOF per control point: u, v, w, θx, θy
%
%           Input : 
%             p,q : polynomial degrees
%              dR : The basis function and up to its second mixed
%                   derivatives
%      GCovariant : = [GCov1 GCov2] the covariant base vectors (2D)
%     dGCovariant : = [dGCov1/dxi dGCov2/deta dGCov1/deta=dGCov2/dxi] the 
%                   derivatives of the covariant base vectors
%         G3Tilde : The not normalized surface normal
%        Dm,Db,Ds : membrane, bending, and shear material matrices
% shearCorrection : Shear correction factor (typically 5/6 for shells)
%
%          Output :
%              Ke : The element stiffness matrix at the Gauss Point
%
% Function layout :
%
% 0. Read input
%
% 1. Compute the covariant metric coefficient tensor
%
% 2. Compute the contravariant base vectors
%
% 3. Compute the local Cartesian basis
%
% 4. Compute the transformation matrix from contravariant to local Cartesian basis
%
% 5. Compute B-operator matrices for membrane, bending, and shear parts
%
% 6. Transform B-operator matrices to local Cartesian system
%
% 7. Compute element stiffness matrices for membrane, bending, and shear parts
%
% 8. Assemble total element stiffness matrix
%
%% Function main body

%% 0. Read input

% The number of DOFs (5 DOF per control point: u, v, w, θx, θy)
numCPsEl = (p + 1)*(q + 1);
numDOFsEl = 5*numCPsEl;

% Initialize auxiliary arrays
BOperatorMatrixMembraneCurvilinear = zeros(3, numDOFsEl);
BOperatorMatrixBendingCurvilinear = zeros(3, numDOFsEl);
BOperatorMatrixShearCurvilinear = zeros(2, numDOFsEl);

% Normalize surface normal
G3 = G3Tilde / norm(G3Tilde);

%% 1. Compute the covariant metric coefficient tensor
GabCov = GCovariant' * GCovariant;

%% 2. Compute the contravariant base vectors
GContravariant = (GabCov \ GCovariant')';

%% 3. Compute the local Cartesian basis
eLC = computeLocalCartesianBasis4BSplineSurface ...
    (GCovariant, GContravariant);

%% 4. Compute the transformation matrix from contravariant to local Cartesian basis
TContra2LC = computeTFromContra2LocalCartesian4VoigtStrainIGAKLShell ...
    (eLC, GContravariant);

%% 5. Compute B-operator matrices for membrane, bending, and shear parts

% Loop over all the degrees of freedom
for iDOFs = 1:numDOFsEl
    % Compute local node number k and dof direction dir
    k = ceil(iDOFs/5);
    dir = iDOFs - 5*(k - 1);
    
    %% 5.1 Membrane strains (in-plane strains)
    % εξξ = ∂u/∂ξ·G1, εηη = ∂v/∂η·G2, γξη = (∂u/∂η·G2 + ∂v/∂ξ·G1)
    
    if dir == 1 % u-displacement
        BOperatorMatrixMembraneCurvilinear(1, iDOFs) = dR(k, 2) * GCovariant(1, 1);  % ∂u/∂ξ·G1x
        BOperatorMatrixMembraneCurvilinear(3, iDOFs) = 0.5 * dR(k, 4) * GCovariant(1, 2); % 0.5*∂u/∂η·G2x
    elseif dir == 2 % v-displacement
        BOperatorMatrixMembraneCurvilinear(2, iDOFs) = dR(k, 4) * GCovariant(2, 2);  % ∂v/∂η·G2y
        BOperatorMatrixMembraneCurvilinear(3, iDOFs) = 0.5 * dR(k, 2) * GCovariant(2, 1); % 0.5*∂v/∂ξ·G1y
    end
    
    %% 5.2 Bending strains (curvature changes)
    % κξξ = -∂θx/∂ξ, κηη = -∂θy/∂η, κξη = -(∂θx/∂η + ∂θy/∂ξ)
    % Note: In HSDT, curvatures are directly related to rotation gradients
    
    if dir == 4 % θx-rotation contributes to bending
        BOperatorMatrixBendingCurvilinear(1, iDOFs) = -dR(k, 2); % -∂θx/∂ξ  
        BOperatorMatrixBendingCurvilinear(3, iDOFs) = -0.5 * dR(k, 4); % -0.5*∂θx/∂η
        
    elseif dir == 5 % θy-rotation contributes to bending
        BOperatorMatrixBendingCurvilinear(2, iDOFs) = -dR(k, 4); % -∂θy/∂η
        BOperatorMatrixBendingCurvilinear(3, iDOFs) = -0.5 * dR(k, 2); % -0.5*∂θy/∂ξ
    end
    
    %% 5.3 Shear strains (transverse shear)
    % γξz = ∂w/∂ξ + θx, γηz = ∂w/∂η + θy
    % Note: In parametric coordinates, shear relates w-gradients to rotations
    
    if dir == 3 % w-displacement contributes to shear
        BOperatorMatrixShearCurvilinear(1, iDOFs) = dR(k, 2); % ∂w/∂ξ
        BOperatorMatrixShearCurvilinear(2, iDOFs) = dR(k, 4); % ∂w/∂η
        
    elseif dir == 4 % θx-rotation contributes to shear  
        BOperatorMatrixShearCurvilinear(1, iDOFs) = dR(k, 1); % θx
        
    elseif dir == 5 % θy-rotation contributes to shear
        BOperatorMatrixShearCurvilinear(2, iDOFs) = dR(k, 1); % θy
    end
end

%% 6. Transform B-operator matrices to local Cartesian system

% For the membrane part (3x3 transformation)
BOperatorMatrixMembraneCartesian = TContra2LC * BOperatorMatrixMembraneCurvilinear;

% For the bending part (3x3 transformation)
BOperatorMatrixBendingCartesian = TContra2LC * BOperatorMatrixBendingCurvilinear;

% For the shear part (2x2 transformation) 
TShear = TContra2LC(1:2, 1:2); % Extract 2x2 submatrix for shear
BOperatorMatrixShearCartesian = TShear * BOperatorMatrixShearCurvilinear;

%% 7. Compute element stiffness matrices for membrane, bending, and shear parts

% Membrane stiffness
Kem = BOperatorMatrixMembraneCartesian' * Dm * BOperatorMatrixMembraneCartesian;

% Bending stiffness  
Keb = BOperatorMatrixBendingCartesian' * Db * BOperatorMatrixBendingCartesian;

% Shear stiffness (with shear correction factor)
Kes = BOperatorMatrixShearCartesian' * (shearCorrection * Ds) * BOperatorMatrixShearCartesian;

%% 8. Assemble total element stiffness matrix
Ke = Kem + Keb + Kes;

end