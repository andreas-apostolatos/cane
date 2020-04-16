function [KTangent, resVct, BSplinePatch, propCoupling, minElASize] = ...
    computeIGAVMSStabMtxAndVct4BossakTINewtonNLinear4StokesE2D ...
    (constMtx, tanMtxLoad, up, upSaved, upDot, upDotSaved, BSplinePatch, ...
    connections, propCoupling, loadFactor, noPatch, noTimeStep, iNLinearIter, ...
    noWeakDBCCnd, t, propFldDynamics, isReferenceUpdated, tab, propInt)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the system matrix and the right-hand side vector corresponding to
% the application of the Bossak time integration scheme to the linear
% isogeometric Stokes problem in 2D which is stabilized using the
% Variational Multiscale Stabilization (VMS) method.
%
%              Input :
%           constMtx : Constant part of the stiffness matrix (dummy)
%         tanMtxLoad : Tangent matrix of follower load (dummy)
%                 up : The discrete solution vector of the previous 
%                      nonlinear iteration step (dummy variable for this 
%                      function)
%            upSaved : The discrete solution vector of the previous time 
%                      step
%              upDot : The rate of the discrete solution vector of the
%                      previous  nonlinear iteration step (dummy variable 
%                      for function)
%         upDotSaved : The rate of the discrete solution vector of the
%                      previous time step
%       BSplinePatch : Polynomial degrees and knot vectors of the 
%                      underlying patch to the isogeometric system
%        connections : Array containing information on the patch coupling 
%                      (dummy)
%       propCoupling : Structure containing information on the patch
%                      coupling (dummy)
%         loadFactor : Load factor in case of load stepping
%            noPatch : dummy
%         noTimeStep : Time step ID (dummy)
%       iNLinearIter : Index of nonlinear iteration step
%       noWeakDBCCnd : dummy
%                  t : Time instance
%    propFldDynamics : Structure containing information on the time
%                      integration regarding the fluid dynamics,
%                              .method : The time integration method
%                           .alphaBeta : (parameter for the Bossak scheme)
%                               .gamma : (parameter for the Bossak scheme)
%                              .TStart : Start time of the simulation
%                                .TEnd : End time of the simulation
%                         .noTimeSteps : Number of time steps
%                                  .dt : Time step
% isReferenceUpdated : Flag on whether the reference configuration is
%                      updated, e.g. for form-finding analysis (dummy)
%                tab : Tabulation for outputting information to the command
%                      window
%             outMsg : Enables outputting information on the command window
%                      if it is set as 'outputEnabled'
%
%             Output :
%           KTangent : The tangent stiffness matrix corresponding to the
%                      Newton linearization of the Bossak discrete equation
%                      system applied to the transient Stokes equations
%             resVct : The residual vector corresponding to the Newton 
%                      linearization of the Bossak discrete equation system 
%                      applied to the transient Stokes equations
%       BSplinePatch : Dummy output for this function
%       propCoupling : Dummy output for this function
%         minElASize : The minimum element area size in the isogeometric 
%                      mesh
%
% Function layout :
%
% 0. Read input
%
% 1. Choose an integration rule
%
% 2. loop over all elements (knot spans)
%
%    2i. Initialize the element area size
%
%   2ii. Compute the determinant of the Jacobian to the transformation from the NURBS space (xi-eta) to the integration domain [-1,1]x[-1,1] 
%
%  2iii. Create the element freedom table
%
%   2iv. Get the element discrete solution vector of the previous Newton iteration step
%
%    2v. Loop over all Quadrature Points for the integration of the element stiffness
%
%        2v.1. Compute the coordinates of the quadrature points at the NURBS domain via mapping
%
%        2v.2. Issue the quadrature weights the compute the quadrature weight needed in the multiplication
%
%        2v.3. Find the correct spans where xi,eta lie in
%
%        2v.4. Compute the NURBS basis functions and their derivatives (possibly also the base vectors) at the quadrature point
%
%        2v.5. Compute the determinant of the Jacobian (and possibly Hessian) to the transformation from the physical space (x-y) to the NURBS parameter space (xi-eta)
%
%        2v.6. Compute the element area on the Gauss Point and sum up the contribution to it
%
%        2v.7. Compute the stabilization parameters for the momentum and the continuity equations
%
%        2v.8. Compute the element stiffness, tangent stiffness, mass matrices and body force vector at the quadrature point multiplying also via the determinants of the Jacobians to the two transformations and the quadrature weight
%
%        2v.9. Add the contribution from the Gauss Point and assemble to the global system
%
%   2vi. Check if the current minimum element area size is smaller than for the previous element
%
% 3. Compute the system matrix corresponding to the Bossak time integration scheme
%
% 4. Compute the right-hand side vector corresponding to the Bossak time integration scheme
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

% Check if the analysis is transient
isTransient = false;
if ~ischar(propFldDynamics)
    if isfield(propFldDynamics, 'timeDependence')
        if strcmp(propFldDynamics.timeDependence, 'transient')
            isTransient = true;
        end
    end
end

% Get all geometric parameters
p = BSplinePatch.p;
q = BSplinePatch.q;
Xi = BSplinePatch.Xi;
Eta = BSplinePatch.Eta;
CP = BSplinePatch.CP;
isNURBS = BSplinePatch.isNURBS;

% Get the flow parameters of the patch
parameters = BSplinePatch.parameters;

% Get the function handle to the computation of the body forces within the
% patch
computeBodyForces = BSplinePatch.computeBodyForces;

% Get DOF numbering of the patch
DOFNumbering = BSplinePatch.DOFNumbering;

% Get the integration properties of the patch
propInt = BSplinePatch.int;

% Initialize flag for the nature of the body forces
isBodyForceNumeric = 0;
isBodyForceIdenticallyZero = 0;
if isnumeric(computeBodyForces)
    isBodyForceNumeric = 1;
    if norm(computeBodyForces) == 0
        isBodyForceIdenticallyZero = 1;
    end
end

% Check input
numKnots_xi = length(Xi);
numKnots_eta = length(Eta);
numCPs_xi = length(CP(:, 1, 1));
numCPs_eta = length(CP(1, :, 1));
checkInputForBSplineSurface ...
    (p, numKnots_xi, numCPs_xi, q, numKnots_eta, numCPs_eta);

% Total number of degrees of freedom for the 2D vector transport equation
numDOFs = 3*numCPs_xi*numCPs_eta;

% Number of DoFs affected the element under study
numDOFsEl = 3*(p + 1)*(q + 1);

% Initialize global stiffness matrix
KLinear = zeros(numDOFs);

% Initialize the global mass matrix
massMtx = zeros(numDOFs);

% Initialize the global body force vector
if ~isBodyForceNumeric && ~isBodyForceIdenticallyZero
    FBody = zeros(numDOFs,1);
else
    FBody = [];
end

% Get the Neumann boundary conditions
NBC = BSplinePatch.NBC;

% Initialize load vector
BSplinePatch.FGamma = zeros(BSplinePatch.noDOFs, 1);

% Compute a nessecary pre-factor for the Bossak time integration scheme
if isTransient
    preFactor = (1 - propFldDynamics.alphaBeta)/propFldDynamics.gamma/ ...
        propFldDynamics.dt;
end

% Initialize minimum element area size
minElementSize = 1e4;
minElASize = minElementSize;

%% 1. Check for the presence of a non-conservative loading and if there exists 
isConservative = true;
tanMtxLoad = 'undefined';
for iNBC = 1:NBC.noCnd
    if NBC.isFollower(iNBC, 1)
        isConservative = false;
        tanMtxLoad = zeros(BSplinePatch.noDOFs);
        BSplinePatch.FNonConservative = zeros(BSplinePatch.noDOFs, 1);
        break;
    end
end
if ~isConservative
    BSplinePatch.FNonConservative = zeros(BSplinePatch.noDOFs, 1);
end

%% 2. Compute the load vector corresponding to boundary applied fluxes
for iNBC = 1:NBC.noCnd
    %% 2i. Initialize the load vector for the current condition
    FGamma = zeros(BSplinePatch.noDOFs,1);

    %% 2ii. Get the function handle for the load vector computation
    funcHandle = str2func(NBC.computeLoadVct{iNBC});

    %% 2iii. Compute the load vector and the tangent matrix resulting from the application of follower loads
    if ~(propFldDynamics.isStaticStep && NBC.isTimeDependent(iNBC, 1))
        [FGamma, tanMtxLoadPatch] = ...
            funcHandle ...
            (FGamma, BSplinePatch, NBC.xiLoadExtension{iNBC}, ...
            NBC.etaLoadExtension{iNBC}, NBC.loadAmplitude{iNBC}, ...
            NBC.loadDirection{iNBC}, NBC.isFollower(iNBC,1), t, ...
            propInt, '');
        if NBC.isFollower(iNBC, 1)
            tanMtxLoad = tanMtxLoad + tanMtxLoadPatch;
        end
    end

    %% 2iv. If the loading is not conservative add the contribution to the non-conservative load vector
    if NBC.isFollower(iNBC, 1)
        BSplinePatch.FNonConservative = BSplinBSplinePatchePatch.FNonConservative + ...
            FGamma;
    end

    %% 2v. Add The compute external load vector into the B-Spline array
    BSplinePatch.FGamma = BSplinePatch.FGamma + FGamma;
end

%% 1. Choose an integration rule

% Get the number of Gauss Points in xi and eta directions
if strcmp(propInt.type,'default')
    numGP_xi = BSplinePatch.p + 1;
    numGP_eta = BSplinePatch.q + 1;
elseif strcmp(propInt.type,'user')
    numGP_xi = propInt.xiNGP;
    numGP_eta = propInt.etaNGP;
end

% Get the Gauss Point coordinates a weights
[xiGP, xiGW] = getGaussPointsAndWeightsOverUnitDomain(numGP_xi);
[etaGP, etaGW] = getGaussPointsAndWeightsOverUnitDomain(numGP_eta);

%% 2. loop over all elements (knot spans)
for j = q + 1:numKnots_eta - q - 1
    for i = p + 1:numKnots_xi - p - 1
        if Xi(i + 1) ~= Xi(i) && Eta(j + 1) ~= Eta(j)
            %% 2i. Initialize the element area size
            
            % Initialize the new one
            minElementSize = 0;
            
            %% 2ii. Compute the determinant of the Jacobian to the transformation from the NURBS space (xi-eta) to the integration domain [-1,1]x[-1,1] 
            
            % This transformation is constant always
            % 
            %         | xi_i+1 - xi_i                    |
            %         | -------------            0       |
            %         |        2                         |
            %  xi,u = |                                  |
            %         |                  eta_j+1 - eta_j |
            %         |        0         --------------- |
            %         |                          2       |
            detJxiu = (Xi(i + 1) - Xi(i))*(Eta(j + 1) - Eta(j))/4;
            
            %% 2iii. Create the element freedom table
            
            % Initialize element freedom table
            EFT = zeros(1, numDOFsEl);
            
            % Initialize counter
            k = 1;
            
            % Relation global-local DoFs
            for cpj = j - q:j
                for cpi = i - p:i
                    EFT(k) = DOFNumbering(cpi, cpj, 1);
                    EFT(k + 1) = DOFNumbering(cpi, cpj, 2);
                    EFT(k + 2) = DOFNumbering(cpi, cpj, 3);
                    
                    % update counter
                    k = k + 3;
                end
            end
            
            %% 2v. Loop over all Quadrature Points for the integration of the element stiffness
            for iGP_xi = 1:length(xiGP)
                for iGP_eta = 1:length(etaGP)
                    %% 2v.1. Compute the coordinates of the quadrature points at the NURBS domain via mapping
                    xi = (Xi(i + 1) + Xi(i) + xiGP(iGP_xi)*(Xi(i + 1) - Xi(i)))/2;
                    eta = (Eta(j + 1) + Eta(j) + etaGP(iGP_eta)*(Eta(j + 1) - Eta(j)))/2;
                
                    %% 2v.2. Issue the quadrature weights the compute the quadrature weight needed in the multiplication
                    GW = xiGW(iGP_xi)*etaGW(iGP_eta);
                
                    %% 2v.3. Find the correct spans where u,v lie in
                    xiSpan = findKnotSpan(xi, Xi, numCPs_xi);
                    etaSpan = findKnotSpan(eta, Eta, numCPs_eta);
                
                    %% 2v.4. Compute the NURBS basis functions and their first and second derivatives
                    dR = computeIGABasisFunctionsAndDerivativesForSurface...
                        (xiSpan, p, xi, Xi, etaSpan, q, eta, Eta, CP, isNURBS, 2);
                        
                    %% 2v.5. Compute the determinant of the Jacobian (and possibly Hessian) to the transformation from the physical space (x-y) to the NURBS parameter space (xi-eta)
                
                    % Initialize Jacobian
                    Jxxi = zeros(2, 2);
                    
                    % Initialize the Hessian matrix
                    Hxxi = zeros(2, 3);
                
                    % initialize counter
                    k = 0;
                
                    % Loop over all the non-zero contributions at the span
                    % under study
                    for c = 0:q
                        for b = 0:p
                            % Update counter
                            k = k + 1;
                        
                            % Compute recursively the entries of the 
                            % Jacobian
                            Jxxi(1, 1) = Jxxi(1, 1) + CP(i - p + b, j - q + c, 1)*dR(k, 2);
                            Jxxi(1, 2) = Jxxi(1, 2) + CP(i - p + b, j - q + c, 2)*dR(k, 2);
                            Jxxi(2, 1) = Jxxi(2, 1) + CP(i - p + b, j - q + c, 1)*dR(k, 4);
                            Jxxi(2, 2) = Jxxi(2, 2) + CP(i - p + b, j - q + c, 2)*dR(k, 4);
                            
                            % Compute recursively the entries of the 
                            % Hessian matrix
                            
                            % d^2 x1/du^2
                            Hxxi(1, 1) = Hxxi(1, 1) + CP(i - p + b, j - q + c, 1)*dR(k, 3);
                            % d^2 x1/dv^2
                            Hxxi(1, 2) = Hxxi(1, 2) + CP(i - p + b, j - q + c, 1)*dR(k, 6);
                            % d^2 x1/dudv
                            Hxxi(1, 3) = Hxxi(1, 3) + CP(i - p + b, j - q + c, 1)*dR(k, 5);
                            % d^2 x2/du^2
                            Hxxi(2, 1) = Hxxi(2, 1) + CP(i - p + b, j - q + c, 2)*dR(k, 3);
                            % d^2 x2/dv^2
                            Hxxi(2, 2) = Hxxi(2, 2) + CP(i - p + b, j - q + c, 2)*dR(k, 6);
                            % d^2 x2/dudv
                            Hxxi(2, 3) = Hxxi(2, 3) + CP(i - p + b, j - q + c, 2)*dR(k, 5);
                        end
                    end
                
                    % Compute the determinant of the Jacobian
                    detJxxi = det(Jxxi);
                    
                    %% 2v.6. Compute the element area on the Gauss Point and sum up the contribution to it
                    minElementSizeOnGP = abs(detJxxi)*abs(detJxiu)*GW;
                    minElementSize = minElementSize + minElementSizeOnGP;
                    
                    %% 2v.7. Compute the stabilization parameters for the momentum and the continuity equations
                    
                    % Stabilization constants
                    CI = 4;
                    Ct = 4;
                    
                    % Stabilization for the momentum equation
                    GMatrix = Jxxi'*Jxxi;
                    if isTransient
                        tauM = (Ct/propFldDynamics.dt^2 + ...
                            CI*parameters.nue^2*norm(GMatrix, 'fro')^2)^(-1/2);
                    else
                        tauM = ...
                            (CI*parameters.nue^2*norm(GMatrix, 'fro')^2)^(-1/2);
                    end

                    % Stabilization for the continuity equation
                    GVector = GMatrix(1, :) + GMatrix(2, :);
                    tauC = (tauM*(GVector*GVector'))^(-1);
                    
                    %% 2v.8. Compute the element stiffness, tangent stiffness, mass matrices and body force vector
                    [KLinearElOnGP, massMtxElOnGP, FBodyElOnGP] =...
                        computeIGAVMSStabElMtxAndVctNLinear4StokesE2D ...
                        (dR, i, p, xi, Xi, j, q, eta, Eta, CP, parameters, ...
                        computeBodyForces, tauM, tauC, Jxxi, Hxxi);

                    %% 2v.9. Add the contribution from the Gauss Point and assemble to the global system
                    
                    % Global linear stiffness matrix of the system
                    KLinear(EFT, EFT) = KLinear(EFT, EFT) + ...
                        KLinearElOnGP*minElementSizeOnGP;
                    
                    % Global mass matrix of the system
                    massMtx(EFT, EFT) = massMtx(EFT, EFT) + ...
                        massMtxElOnGP*minElementSizeOnGP;
                    
                    % For the body force vector of the system
                    if ~isBodyForceNumeric && ~isBodyForceIdenticallyZero
                        FBody(EFT) = FBody(EFT) + FBodyElOnGP*minElementSizeOnGP;
                    end
                end
            end
    
            %% 2vi. Check if the current minimum element area size is smaller than for the previous element
            if minElementSize <= minElASize
                minElASize = minElementSize;
            end
        end
    end
end

%% 3. Compute the system matrix corresponding to the Bossak time integration scheme
if isTransient
    KTangent = preFactor*massMtx + KLinear;
else
    KTangent = KLinear;
end

%% 4. Compute the right-hand side vector corresponding to the Bossak time integration scheme
if isTransient
    if ~isBodyForceNumeric && ~isBodyForceIdenticallyZero
        resVct = (preFactor*massMtx + KLinear)*up - (FBody + BSplinePatch.FGamma) - ...
            ((1 - propFldDynamics.alphaBeta)/propFldDynamics.gamma/ ...
            propFldDynamics.dt)*massMtx*upSaved - ...
            ((1 - propFldDynamics.alphaBeta)/propFldDynamics.gamma - 1)*massMtx*upDotSaved;
    else
        resVct = (preFactor*massMtx + KLinear)*up - (BSplinePatch.FGamma) - ...
            ((1 - propFldDynamics.alphaBeta)/propFldDynamics.gamma/ ...
            propFldDynamics.dt)*massMtx*upSaved - ...
            ((1 - propFldDynamics.alphaBeta)/propFldDynamics.gamma - 1)*massMtx*upDotSaved;
    end
else
    if ~isBodyForceNumeric && ~isBodyForceIdenticallyZero
        resVct = BSplinePatch.FGamma + FBody;
    else
        resVct = BSplinePatch.FGamma;
    end
end 

%% 5. Check output
if isBSplinePatchCell
    BSplinePatch = {BSplinePatch};
end

end