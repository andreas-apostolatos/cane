function [tangMtxWeakDBCNitsche, resVctWeakDBCNitsche, BSplinePatch] = ...
    computeWeakDBCTangMtxResVctNitscheIGAMembrane ...
    (BSplinePatch, dHat, connections, numDOFs, propCoupling, ...
    tanStiffMtx, noPatch, noTimeStep, noNonlinearIteration, ...
    noWeakDBCCond, thickness, t, propTransientAnalysis, tab, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the complement of the stiffness matrix corresponding to an
% isogeometric membrane which accounts for the application of the Dirichlet
% boundary conditions weakly using the Nitsche method. Depending on
% whethear stabilization is automatically estimated or not, stabilization
% terms are added into the variational problem by solving an eigenvalue
% problem at the first nonlinear iteration step only.
%
%                Input :
%         BSplinePatch : The B-Spline patch array containing:
%                       .p,q : The polynomial orders of the B-Spline 
%                              surface in both parametric directions
%                    .Xi,Eta : The knot vectors in both parametric
%                              directions
%                        .CP : The set of control points and weights
%                   .isNURBS : Flag on whether the basis is a NURBS or a
%                              B-Spline
%                .parameters : Technical parameters for the structure
%                   .homDOFs : The global numbering of the DOFs where 
%                              homogeneous Dirichlet boundary conditions 
%                              are applied
%                 .inhomDOFs : The global numbering of the DOFs where
%                              inhomogeneous Dirichlet boundary conditions
%                              are applied
%           .valuesInhomDOFs : The values on the DOFs corresponding to the
%                              application of inhomogeneous Dirichlet
%                              boundary conditions
%                       .NBC : Structure containing information on the 
%                              application of the Neumann boundary 
%                              conditions
%                                        .noCnd : Number of Neumann 
%                                                 boundary conditions
%                              .xiLoadExtension : Cell array {.noCnd} 
%                                                 containing the load 
%                                                 extensions in the xi-
%                                                 direction
%                             .etaLoadExtension : Cell array {.noCnd} 
%                                                 containing the load 
%                                                 extensions in the eta-
%                                                 direction
%                                .loadAmplitude : Array (1,.noCnd) 
%                                                 containing the load 
%                                                 amplitudes
%                                .loadDirection : Array (1,.noCnd) 
%                                                 containing the load 
%                                                 directions
%                               .computeLoadVct : Cell array {.noCnd} 
%                                                 containing the function 
%                                                 name for the computation 
%                                                 of the load vector
%                               .isConservative : Array (1,.noCnd) of flags 
%                                                 indicating whether the 
%                                                 load is conservative or 
%                                                 not
%                   .weakDBC : Structure containing information on the
%                              application of weak boundary conditions:
%                                        .alpha : The penalty parameter
%                                  .xiExtension : xi-extension of the
%                                                 boundary where the weak
%                                                 boundary conditions are
%                                                 applied
%                                 .etaExtension : eta-extension of the
%                                                 boundary where the weak
%                                                 boundary conditions are
%                                                 applied
%                       dHat : The displacement solution vector from the 
%                              previous nonlinear iteration step
%                connections : Dummy variable for this function
%                    numDOFs : Dummy variable for this function
%               propCoupling : Dummy variable for this function
%                tanStiffMtx : The tangent stiffness matrix at the current 
%                              nonlinear iteration step
%                    noPatch : Patch number in the multipatch array
%                 noTimeStep : Number of time step
%       noNonlinearIteration : Number of the nonlinear iteration
%              noWeakDBCCond : Number of weak Dirichlet boundary conditions
%                  thickness : The thickness of the membrane
%                          t : The time instance
%      propTransientAnalysis : Structure on the transient analysis :
%                                   .timeDependence : 'true' or 'false'
%                                      .noTimeSteps : Number of time steps
%                        tab : Tabulation related to outputting information 
%                              on the command window
%                     outMsg : Enables outputting information onto the 
%                              command window when chosen as 
%                              'outputEnabled'
%
%                     Output :
%      tangMtxWeakDBCNitsche : The tangent matrix corresponding to the
%                              Nitsche method for the weak application of
%                              boundary conditions
%       resVctWeakDBCNitsche : The residual vector corresponding to the
%                              Nitsche method for the application of weak
%                              boundary conditions
%               BSplinePatch : The updated with the stabilization factors
%                              BSpline patch array
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the conditions
% ->
%    1i. Get the extensions of the Dirichlet boundary where weak boundary conditions using the Nitsche method are to be applied
%
%   1ii. Get the DOFs on the current Dirichlet boundary
%
%  1iii. Get the running and the fixed parameters on the patch boundary where the weak boundary conditions are applied
%
%   1iv. Get the parameter space of the application of the weak boundary conditions using the knot vector information
%
%    1v. Issue Gauss Point coordinates and weights
%
%   1vi. Loop over the elements on the parameter space where the weak boundary conditions are applied
%   ->
%        1vi.1. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
%
%        1vi.2. Loop over the Gauss points
%        ->
%               1vi.2i. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
%
%              1vi.2ii. Compute the NURBS basis functions
%
%             1vi.2iii. Create the element freedom table
%
%              1vi.2iv. Get the displacement vector of the previous iteration step at the Gauss point
%
%               1vi.2v. Compute the matrices containing the basis functions and their derivatives
%
%              1vi.2vi. Compute the covariant base vectors of the reference configuration
%
%             1vi.2vii. Get the motion vector for the given time instance
%
%            1vi.2viii. Compute the residual of the Dirichlet boundary condition
%
%             1vi.2vix. Compute the normal to the boundary vector and transform it into the contavariant basis
%
%              1vi.2vx. Compute the covariant base vectors of the current configuration
%
%              1vi.2xi. Compute the covariant metric coefficients of the current configuration
%
%             1vi.2xii. Compute the contravariant basis
%
%            1vi.2xiii. Compute the local Cartesian basis
%
%             1vi.2xiv. Compute the transformation matrix from the contravariant basis to the local Cartesian one
%
%              1vi.2xv. Compute the transformation matrix from the local Cartesian basis to the covariant one
%
%             1vi.2xvi. Compute the Green-Lagrange strain in the contravariant basis
%
%            1vi.2xvii. Transform the Green-Lagrange strain in the local Cartesian basis
%
%           1vi.2xviii. Compute the prestress values on the local Cartesian coordinate system
%
%             1vi.2xix. Compute the 2nd Piola-kirchhoff stress in the local Cartesian system
%
%              1vi.2xx. Transform the 2nd Piola-kirchhoff stress in the covariant system
%
%             1vi.2xxi. Compute the stress components
%
%            1vi.2xxii. Compute the traction vector
%
%           1vi.2xxiii. Compute the first variation of the Green-Lagrange strain in the contravariant basis with respect to the DOFs
%
%            1vi.2xxiv. Transform the first variation of the Green-Lagrange strain at the local Cartesian basis
%
%             1vi.2xxv. Compute the first variation of the 2nd Piola-Kichhoff stress in the local Cartesian basis
%
%            1vi.2xxvi. Transform the first variation of the 2nd Piola-Kichhoff stress at the covariant basis
%
%           1vi.2xxvii. Compute the first variation of the traction vector
%
%          1vi.2xxviii. Compute the necessary products needed for the second variation of the traction vector
%
%            1vi.2xxix. Compute the determinant of the Jacobian of the transformation from the physical to the parent domain on the GP
%
%             1vi.2xxx. Compute the element length at the GP
%
%            1vi.2xxxi. Compute the element tangent stiffness matrix contribution corresponding to the weak application of the Dirichlet boundary conditions with the Nitsche method and add it to the global matrix
%
%           1vi.2xxxii. Compute the element residual vector contribution corresponding to the weak application of the Dirichlet boundary conditions with the Nitsche method and add it to the global vector
%
%          1vi.2xxxiii. Compute the mass matrix necessary for the stabilization of the variational problem
%
%           1vi.2xxxiv. Computate the second variation of the traction vector if estimation of the stabilization parameter is chosen
%
%            1vi.2xxxv. Compute the Q matrix necessary for the eigenvalue problem related to the estimation of the stabilization parameter
%        <-
%
%  1vii. Estimate the stabilization parameter for the current condition
%
% 1viii. Update the tangent stiffness matrix and the residual vector with the stabilization terms
% <-
%
%% Function main body

%% 0. Read input

% Check the NURBS geometry input
if iscell(BSplinePatch)
    if length(BSplinePatch) > 1
        error('Multipatch NURBS surface is given as input to the computation of the stiffness matrix for a single patch NURBS surface');
    else
        BSplinePatch = BSplinePatch{1};
    end
end

% Check if the time step number is defined and if not assign it to 1
if ischar(noTimeStep)
    noTimeStep = 1;
end

% Tolerance for the determinant of the covariant metric coefficient tensor
scaleTolDet = 1e-2;

% Reassign the analysis arrays
p = BSplinePatch.p;
q = BSplinePatch.q;
Xi = BSplinePatch.Xi;
Eta = BSplinePatch.Eta;
CP = BSplinePatch.CP;
CPd = BSplinePatch.CPd;
isNURBS = BSplinePatch.isNURBS;
parameters = BSplinePatch.parameters;
prestress = parameters.prestress;

% Get the DOF numbering
DOFNumbering = BSplinePatch.DOFNumbering;

% Compute the material matrix
materialMtxVoigt = parameters.E*parameters.t/(1-parameters.nue^2)*...
        [1              parameters.nue 0
         parameters.nue 1              0
         0              0              (1-parameters.nue)/2];

% On the weak Dirichlet boundary conditions
if BSplinePatch.weakDBC.estimationStabilPrm
    if ~isfield(BSplinePatch.weakDBC,'automaticStabilization')
        if ~ischar(propTransientAnalysis) 
            if isfield(propTransientAnalysis,'noTimeSteps')
                BSplinePatch.weakDBC.automaticStabilization = ...
                    zeros(BSplinePatch.weakDBC.noCnd,propTransientAnalysis.noTimeSteps + 1);
            else
                BSplinePatch.weakDBC.automaticStabilization = ...
                    zeros(BSplinePatch.weakDBC.noCnd,1);
            end
        end
    end
end

% Number of Control Points in xi-,eta- directions
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));

% Number of local DOFs
noCPsEl = (p+1)*(q+1);
noDOFsEl = 3*noCPsEl;

% Number of DOFs
numDOFs = BSplinePatch.noDOFs;

% Initialize matrices related to the stabilization of the variational problem
if BSplinePatch.weakDBC.estimationStabilPrm
    massMtx = zeros(numDOFs,numDOFs);
    RMtxG2 = zeros(numDOFs,1);
    if noNonlinearIteration == 1
        QMtx = zeros(numDOFs,numDOFs);
    end
end

% Initialize the element freedom table
EFT = zeros(1,noDOFsEl);

% Initialize auxiliary arrays
RMtx = zeros(3,noDOFsEl);
dRdxiMtx = zeros(3,noDOFsEl);
dRdetaMtx = zeros(3,noDOFsEl);

% Initialize the output arrays
tangMtxWeakDBCNitsche = zeros(numDOFs,numDOFs);
resVctWeakDBCNitsche = zeros(numDOFs,1);

%% 1. Loop over the conditions
for iCnd = 1:BSplinePatch.weakDBC.noCnd
    %% 1i. Get the extensions of the Dirichlet boundary where weak boundary conditions using the Nitsche method are to be applied
    xiExtension = BSplinePatch.weakDBC.xiExtension{iCnd};
    etaExtension = BSplinePatch.weakDBC.etaExtension{iCnd};
    isWeakDBCOverPoint = false;
    
    %% 1ii. Get the DOFs on the current Dirichlet boundary
    weakDBCDOFs = [];
    for direction = 1:3
        weakDBCDOFs = findDofs3D(weakDBCDOFs,xiExtension,etaExtension,direction,CP);
    end
    
    %% 1iii. Get the running and the fixed parameters on the patch boundary where the weak boundary conditions are applied
    if etaExtension(1) == etaExtension(2) && xiExtension(1) ~= xiExtension(2)
        % Coupled region in xi-direction
        couplingRegion = xiExtension;

        % Find the correct spans for the coupled region
        spanStart = findKnotSpan(couplingRegion(1),Xi,nxi);
        spanEnd = findKnotSpan(couplingRegion(2),Xi,nxi) + 1;

        % Corresponding to the coupled region knot span
        weakDBCRgionOnKnotVector = Xi(spanStart:spanEnd);

        % Fixed parameter on the parametric net
        eta = etaExtension(1);

        % Find the span where xiEta it lies in
        etaSpan = findKnotSpan(eta,Eta,neta);

        % Get distances of the corresponding bounding Control Points to the 
        % boundary
        if etaExtension(1) == Eta(1,1)
            distCP = norm(squeeze(CP(1,1,1:3)) - squeeze(CP(end,1,1:3)));
        elseif etaExtension(1) == Eta(1,end)
            distCP = norm(squeeze(CP(1,end,1:3)) - squeeze(CP(end,end,1:3)));
        else
            error('The knot vector Eta and the corresponding boundary where Dirichlet boundary conditions are applied do not match');
        end
        tolDet = scaleTolDet*distCP;
        if distCP == 0
            warning('The edge where weak Dirichlet conditions are applied is collapsed');
        end

        
        % Flag on whether the coupling line is over xi
        isOnXi = true;
    elseif xiExtension(1) == xiExtension(2) && etaExtension(1) ~= etaExtension(2)
        % Coupled region in eta-direction
        couplingRegion = etaExtension;

        % Find the correct spans for the coupled region
        spanStart = findKnotSpan(couplingRegion(1),Eta,neta);   
        spanEnd = findKnotSpan(couplingRegion(2),Eta,neta) + 1;

        % Corresponding to the coupled region knot span
        weakDBCRgionOnKnotVector = Eta(spanStart:spanEnd);

        % Fixed parameter on the parametric net
        xi = xiExtension(1);

        % Find the span where uv it lies in
        xiSpan = findKnotSpan(xi,Xi,nxi);
        
        % Get distances of the corresponding bounding Control Points to the 
        % boundary
        if xiExtension(1) == Xi(1,1)
            distCP = norm(squeeze(CP(1,1,1:3)) - squeeze(CP(1,end,1:3)));
        elseif xiExtension(1) == Xi(1,end)
            distCP = norm(squeeze(CP(end,1,1:3)) - squeeze(CP(end,end,1:3)));
        else
            error('The knot vector Eta and the corresponding boundary where Dirichlet boundary conditions are applied do not match');
        end
       tolDet = scaleTolDet*distCP;
        if distCP == 0
            warning('The edge where weak Dirichlet conditions are applied is collapsed');
        end

        % Flag on whether the coupling line is over eta
        isOnXi = false;
    elseif etaExtension(1) == etaExtension(2) && xiExtension(1) == xiExtension(2)
        % Flag on whether the integration domain is a point
        isWeakDBCOverPoint = true;
        
        % Fixed parameters on the parametric net
        xi = xiExtension(1);
        eta = etaExtension(1);
        
        % Find the correct spans for the point
        xiSpan = findKnotSpan(xi,Xi,nxi);
        etaSpan = findKnotSpan(eta,Eta,neta);
        
        % Corresponding to the coupled region knot span
        weakDBCRgionOnKnotVector = zeros(3,1);
    else
        error('Boundary over which weak Dirichlet boundary conditions are to be imposed wrongly defined');
    end
    
    %% 1iv. Get the parameter space of the application of the weak boundary conditions using the knot vector information
    if ~isWeakDBCOverPoint
        weakDBCRgionOnKnotVector = unique(weakDBCRgionOnKnotVector);
    end

    %% 1v. Issue Gauss Point coordinates and weights
    if ~isWeakDBCOverPoint
        if strcmp(BSplinePatch.weakDBC.int.type,'default')
            if isOnXi
                pDegree = p + 1;
            else
                pDegree = q + 1;
            end
            noGPs = ceil((pDegree + 1)/2);
        elseif strcmp(BSplinePatch.weakDBC.int.type,'user')
            noGPs = BSplinePatch.weakDBC.int.noGPs;
        end
        [GP,GW] = getGaussPointsAndWeightsOverUnitDomain(noGPs);
    else
        noGPs = 1;
    end

    %% 1vi. Loop over the elements on the parameter space where the weak boundary conditions are applied
    for i = 1:length(weakDBCRgionOnKnotVector)-1
        if (weakDBCRgionOnKnotVector(i) ~= weakDBCRgionOnKnotVector(i+1)) || isWeakDBCOverPoint
            %% 1vi.1. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
            if ~isWeakDBCOverPoint
                detJxizeta = (weakDBCRgionOnKnotVector(i+1) - weakDBCRgionOnKnotVector(i))/2;
            end

            %% 1vi.2. Loop over the Gauss points
            for j = 1:noGPs
                %% 1vi.2i. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
                if ~isWeakDBCOverPoint
                    xiEta = ((1-GP(j))*weakDBCRgionOnKnotVector(i) + (1+GP(j))*weakDBCRgionOnKnotVector(i+1))/2;
                end

                %% 1vi.2ii. Compute the NURBS basis functions
                if ~isWeakDBCOverPoint
                    if isOnXi
                        xi = xiEta;
                        xiSpan = findKnotSpan(xi,Xi,nxi);
                    else
                        eta = xiEta;
                        etaSpan = findKnotSpan(eta,Eta,neta);
                    end
                end
                dR = computeIGABasisFunctionsAndDerivativesForSurface...
                        (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,1);

                %% 1vi.2iii. Create the element freedom table
                r = 1;
                for cpj = etaSpan-q:etaSpan
                    for cpi = xiSpan-p:xiSpan
                        EFT(r)   = DOFNumbering(cpi,cpj,1);
                        EFT(r+1) = DOFNumbering(cpi,cpj,2);
                        EFT(r+2) = DOFNumbering(cpi,cpj,3);
                        r = r + 3;
                    end
                end

                %% 1vi.2iv. Get the displacement vector of the previous iteration step at the Gauss point
                dHatEl = dHat(EFT);
                dispVct = computePostprocDisplacementIGAKirchhoffLoveShell...
                    (p,q,dR(:,1),dHatEl);

                %% 1vi.2v. Compute the matrices containing the basis functions and their derivatives
                for iCPs = 1:noCPsEl
                    % R
                    RMtx(1,3*iCPs - 2) = dR(iCPs,1);
                    RMtx(2,3*iCPs - 1) = dR(iCPs,1);
                    RMtx(3,3*iCPs) = dR(iCPs,1);

                    % dR/dxi
                    dRdxiMtx(1,3*iCPs - 2) = dR(iCPs,2);
                    dRdxiMtx(2,3*iCPs - 1) = dR(iCPs,2);
                    dRdxiMtx(3,3*iCPs) = dR(iCPs,2);

                    % dR/deta
                    dRdetaMtx(1,3*iCPs - 2) = dR(iCPs,3);
                    dRdetaMtx(2,3*iCPs - 1) = dR(iCPs,3);
                    dRdetaMtx(3,3*iCPs) = dR(iCPs,3);
                end

                %% 1vi.2vi. Compute the covariant base vectors of the reference configuration
                [A1,A2] = computeBaseVectorsAndDerivativesForBSplineSurface...
                    (xiSpan,p,etaSpan,q,CP,0,dR);
                A3 = cross(A1,A2)/norm(cross(A1,A2));
                
                %% 1vi.2vii. Get the motion vector for the given time instance
                if isfield(BSplinePatch.weakDBC,'imposedMotion')
                    if isa(BSplinePatch.weakDBC.imposedMotion{iCnd},'function_handle')
                        % Find the Cartesian coordinates of the surface
                        % parameters
                        X = computeCartesianCoordinatesOfAPointOnBSplineSurface...
                            (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,dR(:,1));
                        
                        % Get the motion at the given Cartesian location
                        % and for the given time instance
                        g2 = BSplinePatch.weakDBC.imposedMotion{iCnd}(X(1,1),X(2,1),X(3,1),t);
                    else
                        error('Specify the variable BSplinePatch.weakDBC.imposedMotion as function handle');
                    end
                else
                    g2 = zeros(3,1);
                end
                
                %% 1vi.2viii. Compute the residual of the Dirichlet boundary condition
                resWeakDBC = dispVct - g2;

                %% 1vi.2ix. Compute the normal to the boundary vector and transform it into the contavariant basis
                if ~isWeakDBCOverPoint
                    [uGC,~] = computeNormalAndTangentVectorsToBSplineBoundary...
                        (xi,Xi,eta,Eta,A1,A2,A3,isOnXi);
                else
                    isOnXi = true;
                    [uGCXi,~] = computeNormalAndTangentVectorsToBSplineBoundary...
                        (xi,Xi,eta,Eta,A1,A2,A3,isOnXi);
                    isOnXi = false;
                    [uGCEta,~] = computeNormalAndTangentVectorsToBSplineBoundary...
                        (xi,Xi,eta,Eta,A1,A2,A3,isOnXi);
%                     uGC = (uGCXi + uGCEta)/2;
                    if i == 1
                        uGC = uGCXi;
                    elseif i == 2
                        uGC = uGCEta;
                    end
                    
%                     XCartesian = computeCartesianCoordinatesOfAPointOnBSplineSurface...
%                         (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,dR(:,1));
%                     color = 'black';
%                     lineWidth = 1;
%                     hold on;
%                     plot_vectorArrow(XCartesian,XCartesian + uGC,color,lineWidth);
                end
                uContravariant = [A1'
                                  A2']*uGC;
                                     
                %% 1vi.2x. Compute the covariant base vectors of the current configuration
                [a1,a2] = computeBaseVectorsAndDerivativesForBSplineSurface...
                    (xiSpan,p,etaSpan,q,CPd,0,dR);

                %% 1vi.2xi. Compute the covariant metric coefficients of the current configuration
                aabCovariant = [a1 a2]'*[a1 a2];

                %% 1vi.2xii. Compute the contravariant basis
                AabCovariant = [A1 A2]'*[A1 A2];
                AContravariant = AabCovariant\[A1 A2]';
                AContravariant = AContravariant';

                %% 1vi.2xiii. Compute the local Cartesian basis
                eLC = computeLocalCartesianBasis4BSplineSurface...
                    ([A1 A2],AContravariant);

                %% 1vi.2xiv. Compute the transformation matrix from the contravariant basis to the local Cartesian one
                TFromContraToLC4VoigtStrain = ...
                    computeTFromContra2LocalCartesian4VoigtStrainIGAKLShell...
                    (eLC,AContravariant);

                %% 1vi.2xv. Compute the transformation matrix from the local Cartesian basis to the covariant one
                TFromLCToCov = computeTFromLocalCartesian2CovariantBasis4BSplineSurface...
                    (eLC,AContravariant);

                %% 1vi.2xvi. Compute the Green-Lagrange strain in the contravariant basis
                EpsilonContra = .5*[aabCovariant(1,1) - AabCovariant(1,1)
                                    aabCovariant(2,2) - AabCovariant(2,2)
                                    aabCovariant(1,2) - AabCovariant(1,2)];

                %% 1vi.2xvii. Transform the Green-Lagrange strain in the local Cartesian basis
                EpsilonLC = TFromContraToLC4VoigtStrain*EpsilonContra;

                %% 1vi.2xviii. Compute the prestress values on the local Cartesian coordinate system
                
                % Check if a user defined coordinate system for the 
                % prestresses is chosen
                isPrestressOverDefinedSystem = false;
                if isfield(parameters.prestress,'computeBaseVectors')
                    if ~isfield(parameters.prestress,'computeParametricCoordinates')
                        error('Function handle parameters.prestress.computeParametricCoordinates has to be defined when defining the prestress over a user-defined coordinate system');
                    end
                    isPrestressOverDefinedSystem = true;
                end
                
                % Compute the convective coordinates of the surface
                if isPrestressOverDefinedSystem || isa(parameters.prestress.voigtVector,'function_handle')
                    X = computeCartesianCoordinatesOfAPointOnBSplineSurface...
                        (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,dR(:,1));
                    theta = parameters.prestress.computeParametricCoordinates(X);
                end
                
                % Compute the transformation matrix from the user defined 
                % coordinate system to the local Cartesian coordinate 
                % system if a user defined coordinate system is chosen
                if isPrestressOverDefinedSystem
                    prestressBaseVct = prestress.computeBaseVectors(theta(1,1),theta(2,1));
                    T2LC = computeT2LocalCartesianBasis(prestressBaseVct,eLC);
                else
                    T2LC = [1 0 0
                            0 1 0
                            0 0 1];
                end
                
                % Compute the prestress values
                if isa(parameters.prestress.voigtVector,'function_handle')
                    pTilde = parameters.prestress.voigtVector(theta);
                else
                    pTilde = parameters.prestress.voigtVector;
                end

                % Transform the vector to the local Cartesian space
                % if defined over a user defined coordinate system
                pTilde = T2LC*pTilde;

                %% 1vi.2xix. Compute the 2nd Piola-kirchhoff stress in the local Cartesian system
                NLC = thickness*pTilde +  materialMtxVoigt*EpsilonLC;

                %% 1vi.2xx. Transform the 2nd Piola-kirchhoff stress in the covariant system
                NCovariant = TFromLCToCov*NLC;

                %% 1vi.2xxi. Compute the stress components
                PalphaBeta = [NCovariant(1,1) NCovariant(3,1)
                              NCovariant(3,1) NCovariant(2,1)];

                %% 1vi.2xxii. Compute the traction vector
                tractionVct = [a1 a2]*PalphaBeta*uContravariant;

                %% 1vi.2xxiii. Compute the first variation of the Green-Lagrange strain in the contravariant basis with respect to the DOFs
                dEpsilonContra = [a1(:,1)'*dRdxiMtx
                                  a2(:,1)'*dRdetaMtx
                                  .5*(a2(:,1)'*dRdxiMtx + a1(:,1)'*dRdetaMtx)];

                %% 1vi.2xxiv. Transform the first variation of the Green-Lagrange strain at the local Cartesian basis
                dEpsilonCartesian = TFromContraToLC4VoigtStrain*dEpsilonContra;

                %% 1vi.2xxv. Compute the first variation of the 2nd Piola-Kichhoff stress in the local Cartesian basis
                dNCartesian = materialMtxVoigt*dEpsilonCartesian;

                %% 1vi.2xxvi. Transform the first variation of the 2nd Piola-Kichhoff stress at the covariant basis
                dNCovariant = TFromLCToCov*dNCartesian;

                %% 1vi.2xxvii. Compute the first variation of the traction vector
                dtractionVct = uContravariant(1,1)*a1*dNCovariant(1,:) + ...
                    uContravariant(2,1)*a2*dNCovariant(2,:) + ...
                    (uContravariant(2,1)*a1 + uContravariant(1,1)*a2)*dNCovariant(3,:) + ...
                    (uContravariant(1,1)*NCovariant(1,1) + uContravariant(2,1)*NCovariant(3,1))*dRdxiMtx + ...
                    (uContravariant(2,1)*NCovariant(2,1) + uContravariant(1,1)*NCovariant(3,1))*dRdetaMtx;
                
                %% 1vi.2xxviii. Compute the necessary products needed for the second variation of the traction vector
                PiMtx = [a1*uContravariant(1,1) a2*uContravariant(2,1) a2*uContravariant(1,1) + a1*uContravariant(2,1)]*...
                    TFromLCToCov*materialMtxVoigt*TFromContraToLC4VoigtStrain;
                products = resWeakDBC'*PiMtx;

                %% 1vi.2xxix. Compute the determinant of the Jacobian of the transformation from the physical to the parent domain on the GP
                if ~isWeakDBCOverPoint
                    if isOnXi
                        detJxxi = norm(A1(:,1));
                    else
                        detJxxi = norm(A2(:,1));
                    end
                    if detJxxi < tolDet
                        continue;
                    end
                end

                %% 1vi.2xxx. Compute the element length at the GP
                if ~isWeakDBCOverPoint
                    elementLengthOnGP = detJxxi*detJxizeta*GW(j);
                else
                    elementLengthOnGP = 1;
                end

                %% 1vi.2xxxi. Compute the element tangent stiffness matrix contribution corresponding to the weak application of the Dirichlet boundary conditions with the Nitsche method and add it to the global matrix
                tangMtxWeakDBCNitsche(EFT,EFT) = tangMtxWeakDBCNitsche(EFT,EFT) - ...
                    (dtractionVct'*RMtx + RMtx'*dtractionVct + ...
                    products(1,1)*(dRdxiMtx'*dRdxiMtx) + ...
                    products(1,2)*(dRdetaMtx'*dRdetaMtx) + ...
                    products(1,3)*.5*(dRdxiMtx'*dRdetaMtx + dRdetaMtx'*dRdxiMtx) + ... % [a1 a2]*ddNalphaBeta*uContravariant
                    uContravariant(1,1)*dRdxiMtx'*resWeakDBC*dNCovariant(1,:) + ...
                    uContravariant(2,1)*dRdetaMtx'*resWeakDBC*dNCovariant(2,:) + ...
                    (uContravariant(2,1)*dRdxiMtx' + uContravariant(1,1)*dRdetaMtx')*resWeakDBC*dNCovariant(3,:) + ... % daCovariant_s'*squeeze(dNalphaBeta(rDOFs,:,:))*uContravariant
                    uContravariant(1,1)*dNCovariant(1,:)'*resWeakDBC'*dRdxiMtx + ...
                    uContravariant(2,1)*dNCovariant(2,:)'*resWeakDBC'*dRdetaMtx + ...
                    dNCovariant(3,:)'*resWeakDBC'*(uContravariant(2,1)*dRdxiMtx + uContravariant(1,1)*dRdetaMtx))*elementLengthOnGP; % daCovariant_r'*squeeze(dNalphaBeta(sDOFs,:,:))*uContravariant

                %% 1vi.2xxxii. Compute the element residual vector contribution corresponding to the weak application of the Dirichlet boundary conditions with the Nitsche method and add it to the global vector
                resVctWeakDBCNitsche(EFT,1) = resVctWeakDBCNitsche(EFT,1) - ...
                    (dtractionVct'*resWeakDBC + RMtx'*tractionVct)*elementLengthOnGP;

                %% 1vi.2xxxiii. Compute the mass matrix necessary for the stabilization of the variational problem
                if BSplinePatch.weakDBC.estimationStabilPrm
                    massMtx(EFT,EFT) = massMtx(EFT,EFT) + (RMtx'*RMtx)*elementLengthOnGP;
                    RMtxG2(EFT,1) = RMtxG2(EFT,1) + (RMtx'*g2)*elementLengthOnGP;
                end

                %% 1vi.2xxxiv. Computate the second variation of the traction vector if estimation of the stabilization parameter is chosen
%                 if noNonlinearIteration == 1 && BSplinePatch.weakDBC.estimationStabilPrm == true
%                     % ddtractionVct = (dtractionVct + ...
%                     %                      [dHatEl'*squeeze(ddtractionVct(1,:,:))' 
%                     %                       dHatEl'*squeeze(ddtractionVct(2,:,:))' 
%                     %                       dHatEl'*squeeze(ddtractionVct(3,:,:))'])
%                     ddtractionVct = dtractionVct + ...
%                         [dHatEl'*(PiMtx(1,1)*(dRdxiMtx'*dRdxiMtx) + PiMtx(1,2)*(dRdetaMtx'*dRdetaMtx) + PiMtx(1,3)*.5*(dRdxiMtx'*dRdetaMtx + dRdetaMtx'*dRdxiMtx))
%                          dHatEl'*(PiMtx(2,1)*(dRdxiMtx'*dRdxiMtx) + PiMtx(2,2)*(dRdetaMtx'*dRdetaMtx) + PiMtx(2,3)*.5*(dRdxiMtx'*dRdetaMtx + dRdetaMtx'*dRdxiMtx))
%                          dHatEl'*(PiMtx(3,1)*(dRdxiMtx'*dRdxiMtx) + PiMtx(3,2)*(dRdetaMtx'*dRdetaMtx) + PiMtx(3,3)*.5*(dRdxiMtx'*dRdetaMtx + dRdetaMtx'*dRdxiMtx))] + ...
%                          [dHatEl'*(uContravariant(1,1)*dRdxiMtx(1,:)'*dNCovariant(1,:) + uContravariant(2,1)*dRdetaMtx(1,:)'*dNCovariant(2,:) + (uContravariant(2,1)*dRdxiMtx(1,:)' + uContravariant(1,1)*dRdetaMtx(1,:)')*dNCovariant(3,:))
%                           dHatEl'*(uContravariant(1,1)*dRdxiMtx(2,:)'*dNCovariant(1,:) + uContravariant(2,1)*dRdetaMtx(2,:)'*dNCovariant(2,:) + (uContravariant(2,1)*dRdxiMtx(2,:)' + uContravariant(1,1)*dRdetaMtx(2,:)')*dNCovariant(3,:))
%                           dHatEl'*(uContravariant(1,1)*dRdxiMtx(3,:)'*dNCovariant(1,:) + uContravariant(2,1)*dRdetaMtx(3,:)'*dNCovariant(2,:) + (uContravariant(2,1)*dRdxiMtx(3,:)' + uContravariant(1,1)*dRdetaMtx(3,:)')*dNCovariant(3,:))] + ...
%                           [dHatEl'*(uContravariant(1,1)*dNCovariant(1,:)'*dRdxiMtx(1,:) + uContravariant(2,1)*dNCovariant(2,:)'*dRdetaMtx(1,:) + dNCovariant(3,:)'*(uContravariant(2,1)*dRdxiMtx(1,:) + uContravariant(1,1)*dRdetaMtx(1,:)))
%                            dHatEl'*(uContravariant(1,1)*dNCovariant(1,:)'*dRdxiMtx(2,:) + uContravariant(2,1)*dNCovariant(2,:)'*dRdetaMtx(2,:) + dNCovariant(3,:)'*(uContravariant(2,1)*dRdxiMtx(2,:) + uContravariant(1,1)*dRdetaMtx(2,:)))
%                            dHatEl'*(uContravariant(1,1)*dNCovariant(1,:)'*dRdxiMtx(3,:) + uContravariant(2,1)*dNCovariant(2,:)'*dRdetaMtx(3,:) + dNCovariant(3,:)'*(uContravariant(2,1)*dRdxiMtx(3,:) + uContravariant(1,1)*dRdetaMtx(3,:)))];
%                 end

                %% 1vi.2xxxv. Compute the Q matrix necessary for the eigenvalue problem related to the estimation of the stabilization parameter
                if noNonlinearIteration == 1 && BSplinePatch.weakDBC.estimationStabilPrm
                    
                    % Compute necessary product
                    productsQ = tractionVct'*PiMtx;
                    
%                     QMtx(EFT,EFT) = QMtx(EFT,EFT) + (ddtractionVct'*ddtractionVct)*elementLengthOnGP;

                    QMtx(EFT,EFT) = QMtx(EFT,EFT) + ...
                        (2*(dtractionVct'*dtractionVct) + ...
                        productsQ(1,1)*(dRdxiMtx'*dRdxiMtx) + ...
                        productsQ(1,2)*(dRdetaMtx'*dRdetaMtx) + ...
                        productsQ(1,3)*.5*(dRdxiMtx'*dRdetaMtx + dRdetaMtx'*dRdxiMtx) + ... % [a1 a2]*ddNalphaBeta*uContravariant
                        uContravariant(1,1)*dRdxiMtx'*tractionVct*dNCovariant(1,:) + ...
                        uContravariant(2,1)*dRdetaMtx'*tractionVct*dNCovariant(2,:) + ...
                        (uContravariant(2,1)*dRdxiMtx' + uContravariant(1,1)*dRdetaMtx')*tractionVct*dNCovariant(3,:) + ... % daCovariant_s'*squeeze(dNalphaBeta(rDOFs,:,:))*uContravariant
                        uContravariant(1,1)*dNCovariant(1,:)'*tractionVct'*dRdxiMtx + ...
                        uContravariant(2,1)*dNCovariant(2,:)'*tractionVct'*dRdetaMtx + ...
                        dNCovariant(3,:)'*tractionVct'*(uContravariant(2,1)*dRdxiMtx + uContravariant(1,1)*dRdetaMtx))*... % daCovariant_r'*squeeze(dNalphaBeta(sDOFs,:,:))*uContravariant
                        elementLengthOnGP;
                end
            end
        end
    end
    
    %% 1vii. Estimate the stabilization parameter for the current condition
    if noNonlinearIteration == 1 && BSplinePatch.weakDBC.estimationStabilPrm
        eigenValues = eig(QMtx(weakDBCDOFs,weakDBCDOFs),full(tanStiffMtx(weakDBCDOFs,weakDBCDOFs)));
        if ~isreal(eigenValues)
            eigenValues = real(eigenValues);
        end
        if ischar(noWeakDBCCond)
            error('Number of weak Dirichlet boundary conditions must be an integer');
        end
        stabilPrm = 4*noWeakDBCCond*max(eigenValues);
        if ~ischar(propTransientAnalysis)
            BSplinePatch.weakDBC.automaticStabilization(iCnd,noTimeStep) = stabilPrm;
        end
        if strcmp(outMsg,'outputEnabled')
            if iCnd == 1
                fprintf(strcat(tab,'Application of weak Dirichlet boundary conditions with the Nitsche method for patch %d : \n'),noPatch);
            end
            fprintf(strcat(tab,'Stabilization for Dirichlet boundary [%d %d]x[%d %d] estimated to %d \n'),...
                xiExtension(1,1),xiExtension(1,2),etaExtension(1,1),etaExtension(1,2),stabilPrm);
            if iCnd == BSplinePatch.weakDBC.noCnd
                fprintf('\n');
            end
        end
    elseif noNonlinearIteration ~= 1 && BSplinePatch.weakDBC.estimationStabilPrm
        if ~ischar(propTransientAnalysis)
            stabilPrm = BSplinePatch.weakDBC.automaticStabilization(iCnd,noTimeStep);
        end
    end
    
    %% 1viii. Update the tangent stiffness matrix and the residual vector with the stabilization terms
    if BSplinePatch.weakDBC.estimationStabilPrm
        tangMtxWeakDBCNitsche = tangMtxWeakDBCNitsche + stabilPrm*massMtx;
        resVctWeakDBCNitsche = resVctWeakDBCNitsche + stabilPrm*massMtx*dHat - stabilPrm*RMtxG2;
    end
end

end
