function [K, F, BSplinePatch, propCoupling, minElArea] = ...
    computeStiffMtxAndLoadVctIGAKirchhoffLoveShellLinear ...
    (KConstant, tanMtxLoad, dHat, dHatSaved, dHatDot, dHatDotSaved, ...
    BSplinePatch, connections, propCoupling, loadFactor, noPatch, ...
    noTimeStep, iNLinearIter, noWeakDBCCnd, t, propStrDynamics, ...
    isReferenceUpdated, tab, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the stiffness matrix and the load vector for the linear 
% Kirchhoff-Love shell formulation and the minimum element area in the mesh 
% for performing a refinement study.
%
%                 Input : 
%             KConstant : Dummy variable for this function
%            tanMtxLoad : Dummy variable for this function
%                  dHat : The displacement field of the previous iteration 
%                        step (dummy variable for this function)
%             dHatSaved : The displacement field of the previous time step 
%                         (dummy variable for this function)
%               dHatDot : The velocity field of the previous iteration step 
%                         (dummy variable for this function)
%          dHatDotSaved : The acceleration field of the previous time step 
%                         (dummy variable for this function)
%
%          BSplinePatch : The B-Spline patch containing the following:
%                          .p,.q : The polynomial degrees
%                       .Xi,.Eta : The knot vectors
%                            .CP : The Control Poin coordinates and weights
%                       .isNURBS : Structure on whether the basis is a 
%                                  BSpline or a NURBS
%                           .NBC : Structure containing information on the
%                                  application of the Neumann boundary
%                                  conditions
%                                   .noCnd : Number of Neumann boundary 
%                                            conditions
%                         .xiLoadExtension : Cell array {.noCnd} containing 
%                                            the load extensions in the xi-
%                                            direction
%                        .etaLoadExtension : Cell array {.noCnd} containing 
%                                            the load extensions in the eta
%                                            -direction
%                           .loadAmplitude : Array (1,.noCnd) containing 
%                                            the load amplitudes
%                           .loadDirection : Array (1,.noCnd) containing 
%                                            the load directions
%                          .computeLoadVct : Cell array {.noCnd} containing 
%                                            the function name for the 
%                                            computation of the load vector
%                          .isConservative : Array (1,.noCnd) of flags 
%                                            indicating whether the load is 
%                                            conservative or not
%                           .int : Structure on the domain integration of 
%                                  the stiffness and body force vector
%                    .parameters : Technical parameters of the shell patch
%                  .DOFNumbering : The numbering of the DOFs within the
%                                  patch itself
%           connections : Dummy variable for this function
%          propCoupling : Structure containing information on the coupling
%                         between patches in a multipatch system (dummy 
%                         variable for this function)
%            loadFactor : Dummy variable for this function
%               noPatch : Dummy variable for this function
%            noTimeStep : Dummy variable for this function
%          iNLinearIter : Dummy variable for this function
%          noWeakDBCCnd : Dummy variable for this function
%                     t : Time instance
%       propStrDynamics : Dummy variable for this function
%    isReferenceUpdated : Dummy variable for this function
%                   tab : Dummy variable for this function
%            loadFactor : Load factor for nonlinear analysis (dummy 
%                         variable for this function)
%                outMsg : Enalbes outputting information onto the cell when
%                         it is chosen as 'outputEnabled' (dummy variable 
%                         for this function)
%
%                Output :
%                     K : master stiffness matrix
%                     F : The externally applied load vector
%        BSplinePatches : Dummy output for this function
%          propCoupling : Dummy output for this function
%             minElArea : The minimum element area in the IGA mesh
%
% Function layout :
%
% 0. Read input
%
% 1. Compute the material matrices
%
% 2. Choose an integration rule
%
% 3. loops over elements
% ->
%    3i. Compute the determinant of the Jacobian to the transformation from the NURBS space (xi-eta) to the integration domain [-1,1]x[-1,1] 
%
%   3ii. Create the Element Freedom Table
%
%  3iii. Initialize the element area
%       
%   3iv. Loop over all the Gauss Points
%   ->
%        3iv.1. Compute the NURBS coordinates xi,eta of the Gauss Point coordinates in the bi-unit interval [-1, 1]
%
%        3iv.2. Compute the NURBS basis function and up to their second derivatives at the Gauss Point
%
%        3iv.3. Compute the covariant base vectors and their first derivatives
%
%        3iv.4. Compute the surface normal (third covariant base vector not normalized)
%
%        3iv.5. Compute the legth of G3Tilde (= area dA)
%
%        3iv.6. Compute the element stiffness matrix at the Gauss point
%
%        3iv.7 Compute the element area on the Gauss Point and add the contribution
%
%        3iv.8. Add the contribution from the Gauss Point
%   <-
%    3v. Find the minimum element area in the mesh
% <-
% 
% 4. Compute the exernally applied load vector
%
% 5. Check output
%
%% Function main body

%% 0. Read input

% Check the NURBS geometry input
isBSplinePatchCell = false;
if iscell(BSplinePatch)
    if length(BSplinePatch) > 1
        error('Multipatch NURBS surface is given as input to the computation of the stiffness matrix for a single patch NURBS surface');
    else
        BSplinePatch = BSplinePatch{1};
    end
    isBSplinePatchCell = true;
end

% Assign back the patch properties
p = BSplinePatch.p;
q = BSplinePatch.q;
Xi = BSplinePatch.Xi;
Eta = BSplinePatch.Eta;
CP = BSplinePatch.CP;
isNURBS = BSplinePatch.isNURBS;
parameters = BSplinePatch.parameters;
int = BSplinePatch.int;
DOFNumbering = BSplinePatch.DOFNumbering;

% Neuman boundary conditions
NBC = BSplinePatch.NBC;

% Number of knots in xi-,eta-direction
numKnots_xi = length(Xi);
numKnots_eta = length(Eta);
numCPs_xi = length(CP(:, 1, 1));
numCPs_eta = length(CP(1, :, 1));

% Check input
checkInputForBSplineSurface(p, numKnots_xi, numCPs_xi, q, numKnots_eta, numCPs_eta);

% Initialize minimum element area in the IGA mesh
tolerance = 1e-4;
if abs(CP(1, 1, 1) - CP(numCPs_xi, 1, 1)) >= tolerance
    minElArea = abs(CP(1, 1, 1) - CP(numCPs_xi, 1, 1));
else
    minElArea = CP(1, 1, 1) - CP(1, numCPs_eta, 1);
end

% Local number of DOFs
numDOFsEl = 3*(p + 1)*(q + 1);

% Number DOFs
numDOFs = 3*numCPs_xi*numCPs_eta;

% Initialize global stiffness matrix
K  = zeros(numDOFs);

%% 1. Compute the material matrices

% Compute the membrane material matrix
Dm = parameters.E*parameters.t/(1 - parameters.nue^2)*...
    [1              parameters.nue 0
	 parameters.nue 1              0
     0              0              (1 - parameters.nue)/2];
                                                 
% Compute the bending material matrix
Db = parameters.E*parameters.t^3/(12*(1 - parameters.nue^2))*...
    [1              parameters.nue 0
     parameters.nue 1              0 
     0              0              (1 - parameters.nue)/2];

%% 2. Choose an integration rule

% Select the integration scheme
if strcmp(int.type, 'default')
    numGP_xi = p + 1;
    numGP_eta = q + 1;
elseif strcmp(int.type, 'user')
    numGP_xi = int.xiNGP;
    numGP_eta = int.etaNGP;
end

% Issue the Gauss Point coordinates and weights
[xiGP, xiGW] = getGaussPointsAndWeightsOverUnitDomain(numGP_xi);
[etaGP, etaGW] = getGaussPointsAndWeightsOverUnitDomain(numGP_eta);

%% 3. loops over elements
for j = q + 1:numKnots_eta - q - 1
    for i = p + 1:numKnots_xi - p - 1
        if Xi(i + 1) ~= Xi(i) && Eta(j + 1) ~= Eta(j)
            %% 3i. Compute the determinant of the Jacobian to the transformation from the NURBS space (xi-eta) to the integration domain [-1,1]x[-1,1] 
            %
            %         | xi_i+1 - xi_i                    |
            %         | -------------            0       |
            %         |        2                         |
            %  xi,u = |                                  |
            %         |                  eta_j+1 - eta_j |
            %         |        0         --------------- |
            %         |                          2       |
            detJxiu = (Xi(i + 1) - Xi(i))*(Eta(j + 1) - Eta(j))/4;
            
            %% 3ii. Create the Element Freedom Table
            
            % Initialize element freedome table
            EFT = zeros(1, numDOFsEl);
            
            % initialize counter
            k = 1;
            
            % relation global-local dof
            for cpj = j - q:j
                for cpi = i - p:i
                    EFT(k) = DOFNumbering(cpi, cpj, 1);
                    EFT(k + 1) = DOFNumbering(cpi, cpj, 2);
                    EFT(k + 2) = DOFNumbering(cpi, cpj, 3);
                    
                    % Update counter
                    k = k + 3;
                end
            end
            
            %% 3iii. Initialize the element area
            elementArea = 0;
            
            %% 3iv. Loop over all the Gauss Points
            for iGP_eta = 1:numGP_eta
                for iGP_xi = 1:numGP_xi
                    %% 3iv.1. Compute the NURBS coordinates u,v of the Gauss Point coordinates in the bi-unit interval [-1, 1]
                    xi = (Xi(i + 1) + Xi(i) + xiGP(iGP_xi)*(Xi(i + 1) - Xi(i)))/2;
                    eta = (Eta(j + 1) + Eta(j) + etaGP(iGP_eta)*(Eta(j + 1) - Eta(j)))/2;
                    
                    %% 3iv.2. Compute the NURBS basis function and up to their second derivatives at the Gauss Point
                    nDrvBasis = 2;
                    dR = computeIGABasisFunctionsAndDerivativesForSurface ...
                        (i, p, xi, Xi, j, q, eta, Eta, CP, isNURBS, nDrvBasis);
                    
                    %% 3iv.3. Compute the covariant base vectors and their first derivatives
                    nDrvBaseVct = 1;
                    [dG1, dG2] = ...
                        computeBaseVectorsAndDerivativesForBSplineSurface ...
                        (i, p, j, q, CP, nDrvBaseVct, dR);
                 
                    %% 3iv.4. Compute the surface normal (third covariant base vector not normalized)
                    G3Tilde = cross(dG1(:, 1), dG2(:, 1));
                    
                    %% 3iv.5. Compute the legth of G3Tilde (= area dA)
                    dA = norm(G3Tilde);
  
                    %% 3iv.6. Compute the element stiffness matrix at the Gauss point
                    KeOnGP = computeElStiffMtxKirchhoffLoveShellLinear ...
                        (p, q, dR, [dG1(:, 1) dG2(:, 1)], [dG1(:, 2) dG2(:, 2) ...
                        dG1(:, 3)], G3Tilde, Dm, Db);
%                     KeOnGP = computeElStiffMtxKirchhoffLoveShellAndreasLinear ...
%                         (p, q, dR, [dG1(:, 1) dG2(:, 1)], [dG1(:, 2) dG2(:, 2) ...
%                         dG1(:, 3)], G3Tilde, Dm, Db);
                    
                    %% 3iv.7 Compute the element area on the Gauss Point and add the contribution
                    elementAreaOnGP = dA*detJxiu*xiGW(iGP_xi)*etaGW(iGP_eta);
                    elementArea = elementArea + elementAreaOnGP;
                    
                    %% 3iv.8. Add the contribution from the Gauss Point
                    K(EFT, EFT) = K(EFT, EFT) + KeOnGP*elementAreaOnGP;
                end
            end
            %% 3v. Find the minimum element area in the mesh
            if elementArea<minElArea
                minElArea = elementArea;
            end
        end
    end
end

%% 4. Compute the exernally applied load vector
if isfield(BSplinePatch, 'FGamma')
    F = BSplinePatch.FGamma;
else
    F = 'undefined';
end

%% 5. Check output
if isBSplinePatchCell
    BSplinePatch = {BSplinePatch};
end

end
