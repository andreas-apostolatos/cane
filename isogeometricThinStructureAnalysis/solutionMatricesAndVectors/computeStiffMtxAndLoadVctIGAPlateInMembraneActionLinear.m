function [stiffMtx, F, BSplinePatch, propCoupling, minElSize] = ...
    computeStiffMtxAndLoadVctIGAPlateInMembraneActionLinear...
    (stiffMtx, tanMtxLoad, dHat, dHatSaved, dHatDot, dHatDotSaved, ...
    BSplinePatch, connections, propCoupling, loadFactor, noPatch, ...
    noTimeStep, noNonlinearIteration, noWeakDBCCnd, t, ...
    propTransientAnalysis, isReferenceUpdated, tab, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the linear master stiffness matrix and the complete right-hand 
% side load vector for an isogeometric plate in membrane plane stress 
% action.
%
%                 Input :
%              stiffMtx : Constant precomputed part of the stiffness matrix
%            tanMtxLoad : Tangent stiffness matrix resulting from the
%                         application of followe loads
%                  dHat : Initial guess for the primary field (dummy input 
%                         for this function)
%             dHatSaved : The discrete solution field of the previous time 
%                         step
%               dHatDot : Initial guess for the time derivative of the 
%                         primary field (dummy input for this function)
%          dHatDotSaved : The time derivative of the discrete solution 
%                         field of the previous time step
%          BSplinePatch : Polynmial orders, knot vectors and control points 
%                         of the B-Spline patch as well as the load vector, 
%                         its technical parameters, Dirichlet boundary 
%                         conditions and integration rule
%                          p,q : polynomial degrees
%                       Xi,Eta : knot vectors
%                           CP : vector containing control point 
%                                coordinates and weights
%                      isNURBS : Flag on whether the geometrical basis is 
%                                NURBS or B-Spline
%                   parameters : Technical and geometrical parameters
%                       FGamma : force vector
%                      homDOFs : Array containing information on 
%                                homogeneous Dirichlet boundary conditions
%                    inhomDOFs : Array containing information on the 
%                                inhomogeneous Dirichlet boundary 
%                                conditions
%              valuesInhomDOFs : Values on the inhomogeneous Dirichlet 
%                                boundary conditions
%                          int : On the numerical integration
%                 DOFNumbering : The global numbering of the DOFs in a 
%                                3-dimensional array
%           connections : Dummy array for this function
%          propCoupling : Dummy array for this function
%            loadFactor : The load factor of the current time step
%            noTimeStep : Number of time step
%  noNonlinearIteration : Number of nonlinear iteration step
%          noWeakDBCCnd : Number of weak Dirichlet boundary conditions
%                      t : The time instance
% propTransientAnalysis : Transient analysis parameters: (dummy variable to 
%                         this function)
%                               .method : Time integration method
%                            .alphaBeta : Bossak parameter
%                                .gamma : Bossak parameter
%                                   .T0 : Start time of the simulation
%                                 .TEnd : End time of the simulation
%                          .noTimeSteps : Number of time steps
%                                   .dt : Time step (numeric or adaptive)
%    isReferenceUpdated : Dummy variable for this function (mainly used in 
%                         form-finding analysis)
%                    tab : Tabulation for the output message in the command
%                          window
%             loadFactor : Dummy variable for this function
%                 outMsg : Whether or not to output message on refinement 
%                          progress 'outputEnabled' : enables output 
%                          information
%
%                 Output :
%               stiffMtx : The master stiffness matrix
%                      F : The complete force vector
%           BSplinePatch : The updated with the stabilization parameters
%                          BSpline patch array when weak enforcement of the
%                          Dirichlet boundary conditions is chosen
%           propCoupling : The updated with the stabilization terms coupling 
%                          properties array
%              minElSize : The minimum element area size
%
% Function layout :
%
% 0. Read input
%
% 1. Choose an integration rule if the stiffness matrix of the patch is not precomputed
%
% 2. loop over all elements (knot spans) if the stiffness matrix of the patch is not precomputed
% ->
%    2i. Compute the determinant of the Jacobian to the transformation from the NURBS space (xi-eta) to the integration domain [-1,1]x[-1,1] 
%
%   2ii. Initialize element area size
%
%  2iii. Loop over all Gauss Points
%  ->
%       2iii.1. Compute the coordinates of the quadrature points at the NURBS domain via mapping
%
%       2iii.2. Compute the quadrature weight as a tensor product
%
%       2iii.3. Find the correct spans where xi,eta lie in
%
%       2iii.4. Compute the IGA basis functions and their first derivatives at the quadrature point
%
%       2iii.5. Compute the determinant of the Jacobian to the transformation from the physical space (x-y) to the NURBS parameter space (xi-eta)
%
%       2iii.6. Compute the element area on the Gauss Point
%
%       2iii.7. Compute the element stiffness matrix at the quadrature point multiplying also via the determinants of the Jacobians to the two transformations and the quadrature weight
%
%       2iii.8. Add the contribution from the Gauss Point and assemble to the global matrices/vectors
%   <-
%
%   2iv. Update the element counter
%
%    2v. Check for the minimum element area size in the mesh
% <-
%
% 3. Compute the exernally applied load vector
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

% Get the patch properties
p = BSplinePatch.p;
q = BSplinePatch.q;
Xi = BSplinePatch.Xi;
Eta = BSplinePatch.Eta;
CP = BSplinePatch.CP;
isNURBS = BSplinePatch.isNURBS;
EFT = BSplinePatch.EFT;
propInt = BSplinePatch.int;

% Compute the number of knots and Control Points in both parametric
% directions
numKnots_xi = length(Xi);
numKnots_eta = length(Eta);
numCPs_xi = length(CP(:, 1, 1));
numCPs_eta = length(CP(1, :, 1));

% Number of global DOFs
numDOFs = 2*numCPs_xi*numCPs_eta;

% Check the precomputed part of the stiffness matrix
isPrecomputed = false;
if isnumeric(stiffMtx)
    numDOFs_precomputed =  length(stiffMtx);
    if numDOFs_precomputed ~= numDOFs
        error('The precomputed part of the stiffness matrix has %d DOFs but the system has %d DOFs', ...
            numDOFs_precomputed, numDOFs);
    else
        isPrecomputed = true;
    end
end
if ~isfield(BSplinePatch, 'isStiffMtxPrecom')
    error('Structure BSplinePatch must define field isStiffMtxPrecom')
else
    if ~islogical(BSplinePatch.isStiffMtxPrecom)
        error('BSplinePatch.isStiffMtxPrecom must be a boolean');
    else
        isStiffMtxPrecom = BSplinePatch.isStiffMtxPrecom;
    end
end
if isStiffMtxPrecom && ~isPrecomputed
    error('Flag BSplinePatch.isStiffMtxPrecom regarding the existence of a precomputed stiffness matrix is true but provided constant part of the stiffness matrix is not provided');
end

% Local number of Control Points
numCPsEl = (p + 1)*(q + 1);

% Number of derivatives for the NURBS basis functions
numDrvs = 1;

% Material matrix D (for plane stress problems)
matMtxVoigt = BSplinePatch.parameters.E/(1 - BSplinePatch.parameters.nue^2)*...
    [1                           BSplinePatch.parameters.nue 0
     BSplinePatch.parameters.nue 1                           0 
     0                           0                           (1 - BSplinePatch.parameters.nue)/2];
 
% Initialize minimum element area size
tol = 1e-4;
if abs(CP(1, 1, 1) - CP(end, 1, 1)) > tol
    minElSize = abs(CP(1, 1, 1) - CP(end, 1, 1));
elseif abs(CP(1, 1, 1) - CP(1, end, 1)) > tol
    minElSize = abs(CP(1, 1, 1) - CP(end, 1, 1));
elseif abs(CP(1, end, 1) - CP(end, end, 1)) > tol
    minElSize = abs(CP(1, end, 1) - CP(end, end, 1));
elseif abs(CP(end, 1, 1) - CP(end, end, 1)) > tol
    minElSize = abs(CP(end, 1, 1) - CP(end, end, 1));
end

% Initialize the element counter
counterEl = 1;

% Initialize output arrays
if ~isStiffMtxPrecom
    stiffMtx = zeros(numDOFs, numDOFs);
end

%% 1. Choose an integration rule if the stiffness matrix of the patch is not precomputed
if ~isStiffMtxPrecom
    if isstruct(propInt)
        if ~isfield(propInt, 'type')
            error('propInt must define a type')
        end
        if strcmp(propInt.type, 'default')
            numGP_xi = p + 1;
            numGP_eta = q + 1;
        elseif strcmp(propInt.type, 'user')
            numGP_xi = propInt.xiNGP;
            numGP_eta = propInt.etaNGP;
        else
            error('Select a valid integration type int.type')
        end
    else
        error('propInt must be a structure');
    end

    % Issue the Gauss Point coordinates and weights
    [GP_xi, GW_xi] = getGaussPointsAndWeightsOverUnitDomain(numGP_xi);
    [GP_eta, GW_eta] = getGaussPointsAndWeightsOverUnitDomain(numGP_eta);
end

%% 2. loop over all elements (knot spans) if the stiffness matrix of the patch is not precomputed
if ~isStiffMtxPrecom
    for iEtaSpan = q + 1:numKnots_eta - q - 1
        for iXiSpan = p + 1:numKnots_xi - p - 1
            % check if we are in a non-zero knot span
            if Xi(iXiSpan + 1) ~= Xi(iXiSpan) && Eta(iEtaSpan + 1) ~= Eta(iEtaSpan)
                %% 2i. Compute the determinant of the Jacobian to the transformation from the NURBS space (xi-eta) to the integration domain [-1,1]x[-1,1] 
                %
                %         | xi_i+1 - xi_i                    |
                %         | -------------            0       |
                %         |        2                         |
                %  xi,u = |                                  |
                %         |                  eta_j+1 - eta_j |
                %         |        0         --------------- |
                %         |                          2       |
                detJxiu = (Xi(iXiSpan + 1) - Xi(iXiSpan))*(Eta(iEtaSpan + 1) - Eta(iEtaSpan))/4;

                %% 2ii. Initialize element area size
                elAreaSize = 0;

                %% 2iii. Loop over all Gauss Points
                for iEta = 1:length(GP_eta)
                    for iXi = 1:length(GP_xi)
                        %% 2iii.1. Compute the coordinates of the quadrature points at the NURBS domain via mapping
                        xi = (Xi(iXiSpan + 1) + Xi(iXiSpan) + GP_xi(iXi)*(Xi(iXiSpan + 1) - Xi(iXiSpan)))/2;
                        eta = (Eta(iEtaSpan + 1) + Eta(iEtaSpan) + GP_eta(iEta)*(Eta(iEtaSpan + 1) - Eta(iEtaSpan)))/2;

                        %% 2iii.2. Compute the quadrature weight as a tensor product
                        GW = GW_xi(iXi)*GW_eta(iEta);

                        %% 2iii.3. Find the correct spans where xi,eta lie in
                        xiSpan = findKnotSpan(xi, Xi, numCPs_xi);
                        etaSpan = findKnotSpan(eta, Eta, numCPs_eta);

                        %% 2iii.4. Compute the IGA basis functions and their first derivatives at the quadrature point
                        dR = computeIGABasisFunctionsAndDerivativesForSurface ...
                            (xiSpan, p, xi, Xi, etaSpan, q, eta, Eta, CP, isNURBS, numDrvs);

                        %% 2iii.5. Compute the determinant of the Jacobian to the transformation from the physical space (x-y) to the NURBS parameter space (xi-eta) and the physical coordinates of the Gauss point

                        % Initialize Jacobian
                        Jxxi = zeros(2, 2);

                        % Initialize Cartesian coordinates
                        x = 0;
                        y = 0;
                        z = 0;

                        % initialize counter
                        k = 0;

                        % Loop over all the non-zero contributions at the span
                        % under study
                        for iJEta = 0:q
                            for iJXi = 0:p
                                % Update counter
                                k = k + 1;

                                % Compute recursively the entries of the Jacobian
                                Jxxi(1, 1) = Jxxi(1, 1) + CP(iXiSpan - p + iJXi, iEtaSpan - q + iJEta, 1)*dR(k, 2);
                                Jxxi(1, 2) = Jxxi(1, 2) + CP(iXiSpan - p + iJXi, iEtaSpan - q + iJEta, 2)*dR(k, 2);
                                Jxxi(2, 1) = Jxxi(2, 1) + CP(iXiSpan - p + iJXi, iEtaSpan - q + iJEta, 1)*dR(k, 3);
                                Jxxi(2, 2) = Jxxi(2, 2) + CP(iXiSpan - p + iJXi, iEtaSpan - q + iJEta, 2)*dR(k, 3);

                                % Compute the Cartesian coordinates of the
                                % Gauss Point
                                x = x + CP(iXiSpan - p + iJXi, iEtaSpan - q + iJEta, 1)*dR(k, 1);
                                y = y + CP(iXiSpan - p + iJXi, iEtaSpan - q + iJEta, 2)*dR(k, 1);
                                z = z + CP(iXiSpan - p + iJXi, iEtaSpan - q + iJEta, 3)*dR(k, 1);
                            end
                        end

                        % Compute the determinant of the Jacobian
                        detJxxi = det(Jxxi);

                        %% 2iii.6. Compute the element area on the Gauss Point
                        elAreaSizeOnGP = abs(detJxxi)*abs(detJxiu)*GW;
                        elAreaSize = elAreaSize + elAreaSizeOnGP;

                        %% 2iii.7. Compute the element stiffness matrix at the quadrature point multiplying also via the determinants of the Jacobians to the two transformations and the quadrature weight
                        Ke = computeElStiffMtxAndLoadVctPlateInMembraneActionLinear ...
                            (numCPsEl, dR(:, 2:3), Jxxi, matMtxVoigt);

                        %% 2iii.8. Add the contribution from the Gauss Point and assemble to the global matrix
                        stiffMtx(EFT(:, counterEl), EFT(:, counterEl)) = ...
                            stiffMtx(EFT(:, counterEl), EFT(:, counterEl)) + Ke*elAreaSizeOnGP;
                    end
                end
                %% 2iv. Update the element counter
                counterEl = counterEl + 1;

                %% 2v. Check for the minimum element area size in the mesh
                if elAreaSize < minElSize
                    minElSize = elAreaSize;
                end
            end
        end
    end
end

%% 3. Compute the exernally applied load vector
if isfield(BSplinePatch, 'FGamma')
    F = BSplinePatch.FGamma;
else
    F = 'undefined';
end

%% 4. Check output
if isBSplinePatchCell
    BSplinePatch = {BSplinePatch};
end

end
