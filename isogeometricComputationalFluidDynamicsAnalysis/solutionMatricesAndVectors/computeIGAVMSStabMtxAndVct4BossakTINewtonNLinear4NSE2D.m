function [KTangent, resVct, BSplinePatch, propCoupling, minElASize] = ...
    computeIGAVMSStabMtxAndVct4BossakTINewtonNLinear4NSE2D ...
    (constMtx, tanMtxLoad, up, upSaved, upDot, upDotSaved, BSplinePatch, ...
    connections, propCoupling, loadFactor, noPatch, noTimeStep, iNLinearIter, ...
    noWeakDBCCnd, t, propFldDynamics, isReferenceUpdated, tab, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the tangent stiffness matrix and the residual vector
% corresponding to the nonlinear equation system occuring from the time
% discretization using the Bossak scheme and the space discretization using
% the isogeometric analysis of the Navier-Stokes equations in 2D. 
%
% Revision on the nonlinear equation system :
%                       
%                          G(u) * u = H                               (1)
%
% Nonlinear equation system (1) can be solved with the Newton method
% provided that the tangent of the system can be computed namely :
%
%               GTangent(u) = G(u) + dG(u)/du * u                     (2)
%
% In equation (2) the term dG(u)/du is obviously a third order tensor. By
% defining the residual of the equation system as :
%
%                   resVct(u) = G(u) * u - H                          (3)
%
% we iterate over all Newton iterations until convergence up to a given
% tolerance has been achieved, namely we solve iteratively the following 
% equation system:
%
%            GTangent(u^(i-1)) * du^i = - resVct(u^(i-1))            (4.1)
%                       u^i = u^(i-1) + du^i                         (4.2)
%
%              Input : 
%           constMtx : Constant part of the stiffness matrix (dummy)
%         tanMtxLoad : Tangent matrix of follower load (dummy)
%                 up : The discrete solution vector of the previous 
%                      nonlinear iteration step
%            upSaved : The discrete solution vector of the previous time 
%                      step
%              upDot : The rate of the discrete solution vector of the
%                      previous  nonlinear iteration step
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
%                      system applied to the transient Navier-Stokes
%                      equations
%             resVct : The residual vector corresponding to the Newton 
%                      linearization of the Bossak discrete equation system 
%                      applied to the transient Navier-Stokes equations
%       BSplinePatch : Updated structure containing information about the
%                      B-Spline patch (dummy)
%       propCoupling : Updated structure containing information on the 
%                      patch coupling (dummy)
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
% 3. Compute the tangent stiffness matrix corresponding to the Bossak time integration scheme
%
% 4. Compute the right-hand side residual vector corresponding to the Bossak time integration scheme
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

% Get DOF numbering of the patch
DOFNumbering = BSplinePatch.DOFNumbering;

% Get the flow parameters of the patch
parameters = BSplinePatch.parameters;

% Get the function handle to the computation of the body forces within the
% patch
computeBodyForces = BSplinePatch.computeBodyForces;

% Get the boundary load vector of the patch
F = BSplinePatch.FGamma;

% Get the integration properties of the patch
propInt = BSplinePatch.int;

% Compute a nessecary pre-factor for the Bossak time integration scheme
if isTransient
    preFactor = (1-propFldDynamics.alphaBeta)/propFldDynamics.gamma/ ...
        propFldDynamics.dt;
end

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
mxi = length(BSplinePatch.Xi);
meta = length(BSplinePatch.Eta);
nxi = length(BSplinePatch.CP(:,1,1));
neta = length(BSplinePatch.CP(1,:,1));
checkInputForBSplineSurface(BSplinePatch.p,mxi,nxi,BSplinePatch.q,meta,neta);

% Total number of degrees of freedom for the 2D vector transport equation
nDOFs = 3*nxi*neta;

% Number of DoFs affected the element under study
nDOFsEl = 3*(BSplinePatch.p+1)*(BSplinePatch.q+1);

% Initialize global stiffness matrix
KLinear = zeros(nDOFs,nDOFs);

% Initialize global tangent stiffness matrix
KNLinear = zeros(nDOFs,nDOFs);

% Initialize the global mass matrix
massMtx = zeros(nDOFs,nDOFs);

% Initialize the global body force vector
if ~isBodyForceNumeric && ~isBodyForceIdenticallyZero
    FBody = zeros(nDOFs,1);
else
    FBody = [];
end

% Initialize minimum element area size
minElementSize = 1e4;
minElASize = minElementSize;

%% 1. Choose an integration rule

% Get the number of Gauss Points in xi and eta directions
if strcmp(propInt.type,'default')
    nGPXi = BSplinePatch.p + 1;
    nGPEta = BSplinePatch.q + 1;
elseif strcmp(propInt.type,'manual')
    nGPXi = propInt.xiNGP;
    nGPEta = propInt.etaNGP;
end

% Get the Gauss Point coordinates a weights
[xiGP,xiGW] = getGaussPointsAndWeightsOverUnitDomain(nGPXi);
[etaGP,etaGW] = getGaussPointsAndWeightsOverUnitDomain(nGPEta);

%% 2. loop over all elements (knot spans)
for j = BSplinePatch.q+1:meta-BSplinePatch.q-1
    for i = BSplinePatch.p+1:mxi-BSplinePatch.p-1
        % check if we are in a non-zero knot span
        if BSplinePatch.Xi(i+1) ~= BSplinePatch.Xi(i) && BSplinePatch.Eta(j+1) ~= BSplinePatch.Eta(j)
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
            detJxiu = (BSplinePatch.Xi(i+1)-BSplinePatch.Xi(i))*(BSplinePatch.Eta(j+1)-BSplinePatch.Eta(j))/4;
            
            %% 2iii. Create the element freedom table
            
            % Initialize element freedom table
            EFT = zeros(1,nDOFsEl);
            
            % Initialize counter
            k = 1;
            
            % Relation global-local DoFs
            for cpj = j-BSplinePatch.q:j
                for cpi = i-BSplinePatch.p:i
                    EFT(k)   = DOFNumbering(cpi,cpj,1);
                    EFT(k+1) = DOFNumbering(cpi,cpj,2);
                    EFT(k+2) = DOFNumbering(cpi,cpj,3);
                    
                    % update counter
                    k = k + 3;
                end
            end
            
            %% 2iv. Get the element discrete solution vector of the previous Newton iteration step
            upe = up(EFT);
            
            %% 2v. Loop over all Quadrature Points for the integration of the element stiffness
            for kxi = 1:length(xiGP)
                for keta =1:length(etaGP)
                    %% 2v.1. Compute the coordinates of the quadrature points at the NURBS domain via mapping
                    xi = ( BSplinePatch.Xi(i+1)+BSplinePatch.Xi(i) + xiGP(keta)*(BSplinePatch.Xi(i+1)-BSplinePatch.Xi(i)) )/2;
                    eta = ( BSplinePatch.Eta(j+1)+BSplinePatch.Eta(j) + etaGP(kxi)*(BSplinePatch.Eta(j+1)-BSplinePatch.Eta(j)) )/2;
                
                    %% 2v.2. Issue the quadrature weights the compute the quadrature weight needed in the multiplication
                    GW = xiGW(keta)*etaGW(kxi);
                
                    %% 2v.3. Find the correct spans where u,v lie in
                    xiSpan = findKnotSpan(xi,BSplinePatch.Xi,nxi);
                    etaSpan = findKnotSpan(eta,BSplinePatch.Eta,neta);
                
                    %% 2v.4. Compute the NURBS basis functions and their first and second derivatives
                    dR = computeIGABasisFunctionsAndDerivativesForSurface...
                        (xiSpan,BSplinePatch.p,xi,BSplinePatch.Xi,etaSpan,BSplinePatch.q,eta,BSplinePatch.Eta,BSplinePatch.CP,BSplinePatch.isNURBS,2);
                        
                    %% 2v.5. Compute the determinant of the Jacobian (and possibly Hessian) to the transformation from the physical space (x-y) to the NURBS parameter space (xi-eta)
                
                    % Initialize Jacobian
                    Jxxi = zeros(2,2);
                    
                    % Initialize the Hessian matrix
                    Hxxi = zeros(2,3);
                
                    % initialize counter
                    k = 0;
                
                    % Loop over all the non-zero contributions at the span
                    % under study
                    for c = 0:BSplinePatch.q
                        for b = 0:BSplinePatch.p
                            % Update counter
                            k = k + 1;
                        
                            % Compute recursively the entries of the 
                            % Jacobian
                            Jxxi(1,1) = Jxxi(1,1) + BSplinePatch.CP(i-BSplinePatch.p+b,j-BSplinePatch.q+c,1)*dR(k,2);
                            Jxxi(1,2) = Jxxi(1,2) + BSplinePatch.CP(i-BSplinePatch.p+b,j-BSplinePatch.q+c,2)*dR(k,2);
                            Jxxi(2,1) = Jxxi(2,1) + BSplinePatch.CP(i-BSplinePatch.p+b,j-BSplinePatch.q+c,1)*dR(k,4);
                            Jxxi(2,2) = Jxxi(2,2) + BSplinePatch.CP(i-BSplinePatch.p+b,j-BSplinePatch.q+c,2)*dR(k,4);
                            
                            % Compute recursively the entries of the 
                            % Hessian matrix
                            
                            % d^2 x1/du^2
                            Hxxi(1,1) = Hxxi(1,1) + BSplinePatch.CP(i-BSplinePatch.p+b,j-BSplinePatch.q+c,1)*dR(k,3);
                            % d^2 x1/dv^2
                            Hxxi(1,2) = Hxxi(1,2) + BSplinePatch.CP(i-BSplinePatch.p+b,j-BSplinePatch.q+c,1)*dR(k,6);
                            % d^2 x1/dudv
                            Hxxi(1,3) = Hxxi(1,3) + BSplinePatch.CP(i-BSplinePatch.p+b,j-BSplinePatch.q+c,1)*dR(k,5);
                            % d^2 x2/du^2
                            Hxxi(2,1) = Hxxi(2,1) + BSplinePatch.CP(i-BSplinePatch.p+b,j-BSplinePatch.q+c,2)*dR(k,3);
                            % d^2 x2/dv^2
                            Hxxi(2,2) = Hxxi(2,2) + BSplinePatch.CP(i-BSplinePatch.p+b,j-BSplinePatch.q+c,2)*dR(k,6);
                            % d^2 x2/dudv
                            Hxxi(2,3) = Hxxi(2,3) + BSplinePatch.CP(i-BSplinePatch.p+b,j-BSplinePatch.q+c,2)*dR(k,5);
                        end
                    end
                
                    % Compute the determinant of the Jacobian
                    detJxxi = det(Jxxi);
                    
                    %% 2v.6. Compute the element area on the Gauss Point and sum up the contribution to it
                    minElementSizeOnGP = abs(detJxxi)*abs(detJxiu)*GW;
                    minElementSize = minElementSize + minElementSizeOnGP;
                    
                    %% 2v.7. Compute the stabilization parameters for the momentum and the continuity equations
                    
                    % Get the solution vector [u v p]' of the previous 
                    % Newton iteration at the Gauss point
                    upElonGP = computeNodalVectorIncompressibleFlow2D(dR(:,1),BSplinePatch.p,BSplinePatch.q,upe);
                    velocityFieldOnGP = upElonGP(1:2);
                    
                    % Stabilization constants
                    CI = 4;
                    Ct = 4;
                    
                    % Stabilization for the momentum equation
                    GMatrix = Jxxi'*Jxxi;
                    if isTransient
                        tauM = (Ct/propFldDynamics.dt^2 + ...
                            velocityFieldOnGP'*GMatrix*velocityFieldOnGP + ...
                            CI*parameters.nue^2*norm(GMatrix, 'fro')^2)^(-1/2);
                    else
                        tauM = (velocityFieldOnGP'*GMatrix*velocityFieldOnGP + ...
                            CI*parameters.nue^2*norm(GMatrix, 'fro')^2)^(-1/2);
                    end

                    % Stabilization for the continuity equation
                    GVector = GMatrix(1,:) + GMatrix(2,:);
                    tauC = (tauM*(GVector*GVector'))^(-1);
                    
                    %% 2v.8. Compute the element stiffness, tangent stiffness, mass matrices and body force vector
                    [KLinearElOnGP,KNLinearElOnGP,massMtxElOnGP,FBodyElOnGP] =...
                        computeIGAVMSStabElMtxAndVctNLinear4NSE2D(dR,i, ...
                        BSplinePatch.p,xi,BSplinePatch.Xi,j,BSplinePatch.q,...
                        eta,BSplinePatch.Eta,BSplinePatch.CP,parameters, ...
                        computeBodyForces,tauM,tauC,upe,upElonGP,Jxxi,Hxxi);

                    %% 2v.9. Add the contribution from the Gauss Point and assemble to the global system
                    
                    % Global linear stiffness matrix of the system
                    KLinear(EFT, EFT) = KLinear(EFT, EFT) + ...
                        KLinearElOnGP*minElementSizeOnGP;
                    
                    % Global nonlinear stiffness matrix of the system
                    KNLinear(EFT, EFT) = KNLinear(EFT, EFT) + ...
                        KNLinearElOnGP*minElementSizeOnGP;
                    
                    % Global mass matrix of the system
                    massMtx(EFT, EFT) = massMtx(EFT, EFT) + ...
                        massMtxElOnGP*minElementSizeOnGP;
                    
                    % For the body force vector of the system
                    if ~isBodyForceNumeric && ~isBodyForceIdenticallyZero
                        FBody(EFT) = FBody(EFT) + ...
                            FBodyElOnGP*minElementSizeOnGP;
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

%% 3. Compute the tangent stiffness matrix corresponding to the Bossak time integration scheme
if isTransient
    KTangent = preFactor*massMtx + KLinear + KNLinear;
else
    KTangent = KLinear + KNLinear;
end

%% 4. Compute the right-hand side residual vector corresponding to the Bossak time integration scheme
if isTransient
    if ~isBodyForceNumeric && ~isBodyForceIdenticallyZero
        resVct = (preFactor*massMtx + KLinear)*up - (FBody + F) - ...
            ((1-propFldDynamics.alphaBeta)/propFldDynamics.gamma/ ...
            propFldDynamics.dt)*massMtx*upSaved - ...
            ((1-propFldDynamics.alphaBeta)/propFldDynamics.gamma - 1)*massMtx*upDotSaved;
    else
        resVct = (preFactor*massMtx + KLinear)*up - F - ...
            ((1-propFldDynamics.alphaBeta)/propFldDynamics.gamma/ ... 
            propFldDynamics.dt)*massMtx*upSaved - ...
            ((1-propFldDynamics.alphaBeta)/propFldDynamics.gamma - 1)*massMtx*upDotSaved;            
    end
else
   if ~isBodyForceNumeric && ~isBodyForceIdenticallyZero
        resVct = F + FBody;
    else
        resVct = F;
    end 
end

%% 5. Check output
if isBSplinePatchCell
    BSplinePatch = {BSplinePatch};
end

end