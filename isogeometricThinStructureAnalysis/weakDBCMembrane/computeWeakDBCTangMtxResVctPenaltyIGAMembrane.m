function [tangMtxWeakDBCPenalty,resVctWeakDBCPenalty,BSplinePatch] = ...
    computeWeakDBCTangMtxResVctPenaltyIGAMembrane...
    (BSplinePatch,dHat,connections,noDOFs,propCoupling,tanStiffMtx,...
    noPatch,noTimeStep,noNonlinearIteration,noWeakDBCCond,...
    thickness,t,propTransientAnalysis,tab,outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the complement of the residual vector corresponding to an
% isogeometric membrane which accounts for the application of the Dirichlet
% boundary conditions weakly using the Penalty method. The function is
% intended for weak application of inhomogeneous Dirichlet boundary
% conditions.
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
%                     noDOFs : Dummy variable for this function
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
%      tangMtxWeakDBCPenalty : Dummy output array for this function
%       resVctWeakDBCPenalty : The residual vector corresponding to the
%                              Penalty method for the application of weak
%                              boundary conditions
%               BSplinePatch : The updated with the stabilization factors
%                              BSpline patch array
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over the conditions
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
%              1vi.2iv. Compute the matrices containing the basis functions
%
%              1vi.2vi. Compute the covariant base vectors of the reference configuration
%
%             1vi.2vii. Get the motion vector for the given time instance
%
%            1vi.2viii. Compute the element length at the GP
%
%             1vii.2ix. Compute the element residual vector contribution corresponding to the weak application of the Dirichlet boundary conditions with the Nitsche method and add it to the global vector
%        <-
%   <-
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

% Reassign the analysis arrays
p = BSplinePatch.p;
q = BSplinePatch.q;
Xi = BSplinePatch.Xi;
Eta = BSplinePatch.Eta;
CP = BSplinePatch.CP;
isNURBS = BSplinePatch.isNURBS;

% Get the DOF numbering
DOFNumbering = BSplinePatch.DOFNumbering;

% Number of Control Points in xi-,eta- directions
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));

% Number of local DOFs
noCPsEl = (p+1)*(q+1);
noDOFsEl = 3*noCPsEl;

% Number of DOFs
noDOFs = BSplinePatch.noDOFs;

% Initialize the element freedom table
EFT = zeros(1,noDOFsEl);

% Initialize auxiliary arrays
RMtx = zeros(3,noDOFsEl);

% Initialize the output arrays
tangMtxWeakDBCPenalty = 'undefined';
resVctWeakDBCPenalty = zeros(noDOFs,1);

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
                
                %% 1vi.2iv. Compute the matrices containing the basis functions
                for iCPs = 1:noCPsEl
                    RMtx(1,3*iCPs - 2) = dR(iCPs,1);
                    RMtx(2,3*iCPs - 1) = dR(iCPs,1);
                    RMtx(3,3*iCPs) = dR(iCPs,1);
                end
                
                %% 1vi.2vi. Compute the covariant base vectors of the reference configuration
                [A1,A2] = computeBaseVectorsAndDerivativesForBSplineSurface...
                    (xiSpan,p,etaSpan,q,CP,0,dR);
                
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
                
                %% 1vi.2vii. Compute the determinant of the Jacobian of the transformation from the physical to the parent domain on the GP
                if ~isWeakDBCOverPoint
                    if isOnXi
                        detJxxi = norm(A1(:,1));
                    else
                        detJxxi = norm(A2(:,1));
                    end
                end
                
                 %% 1vi.2viii. Compute the element length at the GP
                if ~isWeakDBCOverPoint
                    elementLengthOnGP = detJxxi*detJxizeta*GW(j);
                else
                    elementLengthOnGP = 1;
                end
                
                %% 1vii.2ix. Compute the element residual vector contribution corresponding to the weak application of the Dirichlet boundary conditions with the Nitsche method and add it to the global vector
                resVctWeakDBCPenalty(EFT,1) = resVctWeakDBCPenalty(EFT,1) - ...
                    BSplinePatch.weakDBC.alpha*(RMtx'*g2)*elementLengthOnGP;
            end
        end
    end
end

end
