function [KNitscheI, CNitscheI, KNitscheJ, resVctNitscheI, resVctNitscheJ, ...
    propCoupling] = computeDDMTangentNitscheMtcesIGAMembrane ...
    (patchI, patchJ, dHatI, dHatJ, isSameOrientation, connections, ...
    propCoupling, tanStiffMtxI, tanStiffMtxJ, noConnection, noTimeStep, ...
    noNonlinearIteration, propTransientAnalysis, tab, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the tangent stiffness and coupling matrix corresponding to the
% multipatch coupling of membrane multipatches using the Nitsche method as
% well as the matrix which given that it is multiplied with the solution
% vector returns the residual due to the application of the Nitsche method.
% The matrices which are accounting for the coupling between the patches
% are transpose with respect to each other.
%
%                     Input : 
%             patchI,patchJ : The patches which share an interface
%                               .xicoup,.etacoup : The coupling extension
%                                                  boundaries for both of
%                                                  the given patches which
%                                                  share an interface
%               dHat1,dHat2 : The displacement field of the master and the
%                             slave patch from the previous nonlinear 
%                             iteration
%         isSameOrientation : Flag on whether the couling surfaces of the 
%                             shells are oriented in the same direction 
%                             over the coupling interface
%               connections : Define the connection between the patches:
%                            .No : Number of connections
%                     .xiEtaCoup : [patchID1 patchID2 xi12 eta12 xi21 eta21
%                                   ...      ...    ...   ...  ...   ...]
%              propCoupling : Properties of the multipatch coupling
%                             .estimationStabilPrm : Boolean on whether
%                                                    estimation of the
%                                                    stabilization factor
%                                                    is assumed
%                                       gammaTilde : Linear combination
%                                                    factor for the
%                                                    definition of the
%                                                    interface traction
%                                            .intC : On the integration of 
%                                                    the coupling interface
% tanStiffMtxI,tanStiffMtxJ : The tangent stiffness matrices for both
%                             patches (only used when estimation of the 
%                             stabilization is chosen)
%              noConnection : Number of the connection treated
%                noTimeStep : Number of time step
%      noNonlinearIteration : Number of the nonlinear iteration step
%     propTransientAnalysis : Structure on the transient analysis :
%                               .timeDependence : 'true' or 'false'
%                                  .noTimeSteps : Number of time steps
%                       tab : Tabulation for outputting information onto
%                             the command window
%                    outMsg : Enabled outputting information onto the
%                             command window when chosen as 'outputEnabled'
%
%                    Output :
%                   KnI,KnJ : Tangent force coupling contribution to the 
%                             tangent matrices for patch I and patch J
%                       CnI : Tangent force coupling contribution to the 
%                             tangent matrix for patch I (for patch J the 
%                             coupling matrix is CnI')
% FResNitsche1,FResNitsche2 : Residual vectors corresponding to the
%                             additional terms from the Nitsche method
%              propCoupling : The updated with the stabilization terms
%                             coupling properties array
%
% Function layout :
%
% 0. Read input
%
% 1. Get the running and the fixed parameters on the patch interface and the coupling region
%
% 2. Compute the merged knot vector from both patches over the interface
%
% 3. Issue Gauss Point coordinates and weights
%
% 4. Loop over the elements on the coupling surface
% ->
%    4i. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
%
%   4ii. Loop over the Gauss points
%   ->
%        4ii.1. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
%
%        4ii.2. Compute the NURBS basis functions
%
%        4ii.3. Create the element freedom tables
%
%        4ii.4. Get the displacement vectors of the previous iteration step at the Gauss point
%
%        4ii.5. Compute the matrices containing the basis functions and their derivatives
%
%        4ii.6. Compute the covariant base vectors of the reference configuration
%
%        4ii.7. Compute the normal to the boundary vector and transform it into the contavariant basis
%
%        4ii.8. Compute the covariant base vectors of the current configuration
%
%        4ii.9. Compute the covariant metric coefficients of the current configuration
%
%        4ii.10. Compute the contravariant bases
%
%        4ii.11. Compute the local Cartesian bases
%
%        4ii.12. Compute the transformation matrix from the contravariant basis to the local Cartesian one
%
%        4ii.13. Compute the transformation matrix from the local Cartesian bases to the covariant one
%
%        4ii.14. Compute the Green-Lagrange strains in the contravariant bases
%
%        4ii.15. Transform the Green-Lagrange strains in the local Cartesian bases
%
%        4ii.16. Compute the prestress values on the local Cartesian coordinate systems
%
%        4ii.17. Compute the 2nd Piola-kirchhoff stresses in the local Cartesian systems
%
%        4ii.18. Transform the 2nd Piola-kirchhoff stresses in the covariant systems
%
%        4ii.19. Compute the stress components
%
%        4ii.20. Compute the traction vectors
%
%        4ii.21. Compute the first variation of the Green-Lagrange strains in the contravariant bases with respect to the DOFs
%
%        4ii.22. Transform the first variation of the Green-Lagrange strains at the local Cartesian bases
%
%        4ii.23. Compute the first variations of the 2nd Piola-Kichhoff stresses in the local Cartesian bases
%
%        4ii.24. Transform the first variations of the 2nd Piola-Kichhoff stresses at the covariant bases
%
%        4ii.25. Compute the first variation of the traction vectors
%
%        4ii.26. Compute the necessary products needed for the second variations of the traction vectors
%
%        4ii.27. Compute the determinant of the Jacobian of the transformation from the physical to the parent domain on the GP
%
%        4ii.28. Compute the element length at the GP
%
%        4ii.29. Compute the element tangent stiffness matrix contributions on the GP and add them to the global matrices
%
%        4ii.30. Compute the element tangent coupling matrices at the GP
%
%        4ii.31. Compute the element residual vector contributions on the GP and add the to the global vectors
%
%        4ii.32. Compute the mass matrices necessary for the stabilization of the variational problem at the tangent stiffness matrices level
%
%        4ii.33. Compute the mass matrices necessary for the stabilization of the variational problem at the tangent coupling matrices level
%
%        4ii.34. Compute the second variation of the traction vectors needed for the interface eigenvalue problem which needs to be solved for the estimation of the stabilization parameters
%
%        4ii.35. Compute the Q matrix necessary for the eigenvalue problem related to the estimation of the stabilization parameter at the tangent stiffness matrices level
%
%        4ii.36. Compute the Q matrix necessary for the eigenvalue problem related to the estimation of the stabilization parameter at the tangent coupling matrices level
%   <-
% <-
%
% 5. Find the DOFs on the current interface boundary if the interface eigenvalue problem for the estimation of the stabilization factor is to be solved
%
% 6. Estimate the stabilization parameter for the current patch coupling if the interface eigenvalue problem for the estimation of the stabilization factor is to be solved
%
% 7. Stabilize accordingly the tangent and the coupling stiffness matrices if the interface eigenvalue problem for the estimation of the stabilization factor is to be solved
% 
%% Function main body

%% 0. Read input

% On the stabilization of the provblem
if propCoupling.estimationStabilPrm == true
    if ~isfield(propCoupling,'automaticStabilization')
        if ~ischar(propTransientAnalysis)
            if isfield(propTransientAnalysis,'noTimeSteps')
                propCoupling.automaticStabilization = ...
                    zeros(connections.No,propTransientAnalysis.noTimeSteps + 1);
            else
                propCoupling.automaticStabilization = ...
                    zeros(connections.No,1);
            end
        end
    end
end

% Check if the time step number is defined and if not assign it to 1
if ischar(noTimeStep)
    noTimeStep = 1;
end

% Get the number of connections
mI = connections.No;

% For patch I :
% _____________

% Reassign the analysis arrays
pI = patchI.p;
qI = patchI.q;
XiI = patchI.Xi;
EtaI = patchI.Eta;
CPI = patchI.CP;
CPdI = patchI.CPd;
isNURBSI = patchI.isNURBS;
parametersI = patchI.parameters;
thicknessI = parametersI.t;
prestressI = parametersI.prestress;
xicoupI = patchI.xicoup;
etacoupI = patchI.etacoup;

% Compute the material matrix
materialMtxVoigtI = parametersI.E*parametersI.t/(1-parametersI.nue^2)*...
    [1               parametersI.nue 0
     parametersI.nue 1               0
     0               0               (1-parametersI.nue)/2];

% Get the DOF numbering
DOFNumberingI = patchI.DOFNumbering;

% Number of Control Points in xi-,eta- directions
nxiI = length(CPI(:,1,1));
netaI = length(CPI(1,:,1));

% Number of local DOFs
noCPsElI = (pI + 1)*(qI + 1);
noDOFsElI = 3*noCPsElI;

% Number of DOFs
noDOFsI = 3*nxiI*netaI;

% For patch J :
% _____________

pJ = patchJ.p;
qJ = patchJ.q;
XiJ = patchJ.Xi;
EtaJ = patchJ.Eta;
CPJ = patchJ.CP;
CPdJ = patchJ.CPd;
isNURBSJ = patchJ.isNURBS;
parametersJ = patchJ.parameters;
thicknessJ = parametersJ.t;
prestressJ = parametersJ.prestress;
xicoupJ = patchJ.xicoup;
etacoupJ = patchJ.etacoup;

% Compute the material matrix
materialMtxVoigtJ = parametersJ.E*parametersJ.t/(1-parametersJ.nue^2)*...
    [1               parametersJ.nue 0
     parametersJ.nue 1               0
     0               0               (1-parametersJ.nue)/2];
 
% Define a basic tolerance value
tolBasic = 1e-2;

% Get a characteristic length for the interface
isOnXiI = false;
if etacoupI(1) == etacoupI(2)
    isOnXiI = true;
end
if isOnXiI
    distCPStart = norm(squeeze(CPI(1,1,1:3)) - squeeze(CPI(end,1,1:3)));
    distCPEnd = norm(squeeze(CPI(1,end,1:3)) - squeeze(CPI(end,end,1:3)));
else
    distCPStart = norm(squeeze(CPI(1,1,1:3)) - squeeze(CPI(1,end,1:3)));
    distCPEnd = norm(squeeze(CPI(end,1,1:3)) - squeeze(CPI(end,end,1:3)));
end
characteristicLength = max([distCPEnd distCPStart]);
characteristicArea = characteristicLength^2;

% Tolerance for the determinant of the transformation from the parameter to
% the Cartesian space
tolDet = tolBasic*characteristicLength;

% Tolerance for the norm of the surface normal
tolSurfaceNormal = tolBasic*characteristicArea;
 
% Get the DOF numbering
DOFNumberingJ = patchJ.DOFNumbering;

% Number of Control Points in xi-,eta- directions
nxiJ = length(CPJ(:,1,1));
netaJ = length(CPJ(1,:,1));

% Number of local DOFs
noCPsElJ = (pJ+1)*(qJ+1);
noDOFsElJ = 3*noCPsElJ;

% Number of DOFs
noDOFsJ = 3*nxiJ*netaJ;

% Initialize the element freedom tables
EFTI = zeros(1,noDOFsElI);
EFTJ = zeros(1,noDOFsElJ);

% Initialize matrices related to the stabilization of the variational problem
if noNonlinearIteration == 1 && propCoupling.estimationStabilPrm == true
    QKMtxI = zeros(noDOFsI,noDOFsI);
    QKMtxJ = zeros(noDOFsJ,noDOFsJ);
    QCMtxI = zeros(noDOFsI,noDOFsJ);
end
if propCoupling.estimationStabilPrm == true
    massMtxKI = zeros(noDOFsI,noDOFsI);
    massMtxKJ = zeros(noDOFsJ,noDOFsJ);
    massMtxCI = zeros(noDOFsI,noDOFsJ);
end

% Initialize auxiliary arrays
RMtxI = zeros(3,noDOFsElI);
RMtxJ = zeros(3,noDOFsElJ);
dRdxiMtxI = zeros(3,noDOFsElI);
dRdxiMtxJ = zeros(3,noDOFsElJ);
dRdetaMtxI = zeros(3,noDOFsElI);
dRdetaMtxJ = zeros(3,noDOFsElJ);

% Initialize the output arrays
KNitscheI = zeros(noDOFsI,noDOFsI);
KNitscheJ = zeros(noDOFsJ,noDOFsJ);
CNitscheI = zeros(noDOFsI,noDOFsJ);
resVctNitscheI = zeros(noDOFsI,1);
resVctNitscheJ = zeros(noDOFsJ,1);

%% 1. Get the running and the fixed parameters on the patch interface and the coupling region

% For patch I :
% _____________

if etacoupI(1) == etacoupI(2)
    % Coupled region in xi-direction
    couplingRegionI = xicoupI;
    
    % Find the correct spans for the coupled region
    spanStartI = findKnotSpan(couplingRegionI(1),XiI,nxiI);
    spanEndI = findKnotSpan(couplingRegionI(2),XiI,nxiI) + 1;
    
    % Corresponding to the coupled region knot span
    couplingRegionOnKnotVectorI = XiI(spanStartI:spanEndI);
    
    % Fixed parameter on the parametric net
    etaI = etacoupI(1);
    
    % Find the span where xiEta it lies in
    etaSpanI = findKnotSpan(etaI,EtaI,netaI);
    
    % Flag on whether the coupling line is over xi
    isOnXiI = true;
else
    % Coupled region in eta-direction
    couplingRegionI = etacoupI;
    
    % Find the correct spans for the coupled region
    spanStartI = findKnotSpan(couplingRegionI(1),EtaI,netaI);   
    spanEndI = findKnotSpan(couplingRegionI(2),EtaI,netaI) + 1;   
    
    % Corresponding to the coupled region knot span
    couplingRegionOnKnotVectorI = EtaI(spanStartI:spanEndI);
    
    % Fixed parameter on the parametric net
    xiI = xicoupI(1);
    
    % Find the span where uv it lies in
    xiSpanI = findKnotSpan(xiI,XiI,nxiI);
    
    % Flag on whether the coupling line is over eta
    isOnXiI = false;
end

% For patch J :
% _____________

if etacoupJ(1) == etacoupJ(2)
	% Coupled region in xi-direction
    couplingRegionJ = xicoupJ;
    
    % Find the correct spans for the coupled region
    spanStartJ = findKnotSpan(couplingRegionJ(1),XiJ,nxiJ);   
    spanEndJ = findKnotSpan(couplingRegionJ(2),XiJ,nxiJ)+1; 
    
    % Corresponding to the coupled region knot span
    couplingRegionOnKnotVectorJ = XiJ(spanStartJ:spanEndJ);
    

    % Fixed parameter on the parametric net
    etaJ = etacoupJ(1);
    
    % Find the span where xiEta it lies in
    etaSpanJ = findKnotSpan(etaJ,EtaJ,netaJ);
    
    % Flag on whether the coupling line is over xi
    isOnXiJ = true;
else
    % Coupled region in eta-direction
    couplingRegionJ = etacoupJ;
    
    % Find the correct spans for the coupled region
    spanStartJ = findKnotSpan(couplingRegionJ(1),EtaJ,netaJ);   
    spanEndJ = findKnotSpan(couplingRegionJ(2),EtaJ,netaJ)+1;
    
    % Corresponding to the coupled region knot span
    couplingRegionOnKnotVectorJ = EtaJ(spanStartJ:spanEndJ);

    % Fixed parameter on the parametric net
    xiJ = xicoupJ(1);
    
    % Find the span where uv it lies in
    xiSpanJ = findKnotSpan(xiJ,XiJ,nxiJ);
    
    % Flag on whether the coupling line is over eta
    isOnXiJ = false;
end

%% 2. Compute the merged knot vector from both patches over the interface

% Check if the coupling regions from both knot vectors coincide, if not
% transform the coupling knot region on the slave patch to comply with that
% of the master patch
% if couplingRegionOnKnotVector2(1) ~= couplingRegionOnKnotVector1(1) || ...
%         couplingRegionOnKnotVector2(end) ~= couplingRegionOnKnotVector1(end)
%     couplingRegionOnKnotVector2 = transformKnotVct...
%         (couplingRegionOnKnotVector2,[couplingRegionOnKnotVector1(1) couplingRegionOnKnotVector1(end)]);
%     if isOnXi2
%         Xi2 = transformKnotVct(Xi2,[couplingRegionOnKnotVector1(1) couplingRegionOnKnotVector1(end)]);
%     else
%         Eta2 = transformKnotVct(Eta2,[couplingRegionOnKnotVector1(1) couplingRegionOnKnotVector1(end)]);
%     end
% end

% Merge the two knot vectors into one for integration purposes
couplingRegionOnKnotVector = mergesorted(couplingRegionOnKnotVectorI,couplingRegionOnKnotVectorJ);

% Delete double entries
couplingRegionOnKnotVector = unique(couplingRegionOnKnotVector);

%% 3. Issue Gauss Point coordinates and weights
if strcmp(propCoupling.intC.type,'default')
    if isOnXiI
        pDegreeI = pI + 1;
    else
        pDegreeI = qI + 1;
    end
    if isOnXiJ
        pDegreeJ = pJ + 1;
    else
        pDegreeJ = qJ + 1;
    end
    noGPs = ceil((pDegreeI + pDegreeJ + 1)/2);
elseif strcmp(propCoupling.intC.type,'user')
    noGPs = propCoupling.intC.noGPs;
end
[GP,GW] = getGaussPointsAndWeightsOverUnitDomain(noGPs);

%% 4. Loop over the elements on the coupling surface
for i = 1:length(couplingRegionOnKnotVector)-1
    if couplingRegionOnKnotVector(i) ~= couplingRegionOnKnotVector(i+1)
        %% 4i. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
        detJxizeta = (couplingRegionOnKnotVector(i+1)-couplingRegionOnKnotVector(i))/2;

        %% 4ii. Loop over the Gauss points
        for j = 1:noGPs
            %% 4ii.1. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
            xiEta = ((1-GP(j))*couplingRegionOnKnotVector(i)+(1+GP(j))*couplingRegionOnKnotVector(i+1))/2;

            %% 4ii.2. Compute the NURBS basis functions
            
            % For patch I :
            % _____________
            
            if isOnXiI
                xiI = xiEta;
                xiSpanI = findKnotSpan(xiI,XiI,nxiI);
            else
                etaI = xiEta;
                etaSpanI = findKnotSpan(etaI,EtaI,netaI);
            end
            dRI = computeIGABasisFunctionsAndDerivativesForSurface...
                (xiSpanI,pI,xiI,XiI,etaSpanI,qI,etaI,EtaI,CPI,isNURBSI,1);
            
            % For patch J :
            % _____________
            
            if isOnXiJ
                xiJ = xiEta;
                if ~isSameOrientation
                    xiJ = XiJ(length(XiJ)) - xiJ;
                end
                xiSpanJ = findKnotSpan(xiJ,XiJ,nxiJ);
            else
                etaJ = xiEta;
                if ~isSameOrientation
                    etaJ = EtaJ(length(EtaJ)) - etaJ;
                end
                etaSpanJ = findKnotSpan(etaJ,EtaJ,netaJ);
            end
            dRJ = computeIGABasisFunctionsAndDerivativesForSurface...
                (xiSpanJ,pJ,xiJ,XiJ,etaSpanJ,qJ,etaJ,EtaJ,CPJ,isNURBSJ,1);
            
            %% Make a check
%             XI = computeCartesianCoordinatesOfAPointOnBSplineSurface...
%                 (xiSpanI,pI,xiI,XiI,etaSpanI,qI,etaI,EtaI,CPI,dRI(:,1));
%             XJ = computeCartesianCoordinatesOfAPointOnBSplineSurface...
%                 (xiSpanJ,pJ,xiJ,XiJ,etaSpanJ,qJ,etaJ,EtaJ,CPJ,dRJ(:,1));
%             if norm(XI - XJ) > 1e-10;
%                 warning('The integration points for the multipatches have different image in the physical space ||diff|| = %d',norm(XI - XJ));
%             end
            
            %% 4ii.3. Create the element freedom tables
            
            % For patch I :
            % _____________
            
            % Initialize of the counter
            rI = 1;

            % Relation global-local DoFs
            for cpj = etaSpanI-qI:etaSpanI
                for cpi = xiSpanI-pI:xiSpanI
                    EFTI(rI)   = DOFNumberingI(cpi,cpj,1);
                    EFTI(rI+1) = DOFNumberingI(cpi,cpj,2);
                    EFTI(rI+2) = DOFNumberingI(cpi,cpj,3);

                    % update counter
                    rI = rI + 3;
                end
            end
            
            % For patch J :
            % _____________
                        
            % Initialize of the counter
            rJ = 1;

            % Relation global-local DoFs
            for cpj = etaSpanJ-qJ:etaSpanJ
                for cpi = xiSpanJ-pJ:xiSpanJ
                    EFTJ(rJ)   = DOFNumberingJ(cpi,cpj,1);
                    EFTJ(rJ+1) = DOFNumberingJ(cpi,cpj,2);
                    EFTJ(rJ+2) = DOFNumberingJ(cpi,cpj,3);

                    % update counter
                    rJ = rJ + 3;
                end
            end
            
            %% 4ii.4. Get the displacement vectors of the previous iteration step at the Gauss point
            
            % For patch I :
            % _____________
            
            dHatElI = dHatI(EFTI);
            dispVctI = computePostprocDisplacementIGAKirchhoffLoveShell...
                (pI,qI,dRI(:,1),dHatElI);
            
            % For patch J :
            % _____________
            
            dHatElJ = dHatJ(EFTJ);
            dispVctJ = computePostprocDisplacementIGAKirchhoffLoveShell...
                (pJ,qJ,dRJ(:,1),dHatElJ);
            
            %% 4ii.5. Compute the matrices containing the basis functions and their derivatives
            
            % For patch I :
            % _____________
            
            for iCPs = 1:noCPsElI
                % R
                RMtxI(1,3*iCPs - 2) = dRI(iCPs,1);
                RMtxI(2,3*iCPs - 1) = dRI(iCPs,1);
                RMtxI(3,3*iCPs) = dRI(iCPs,1);

                % dR/dxi
                dRdxiMtxI(1,3*iCPs - 2) = dRI(iCPs,2);
                dRdxiMtxI(2,3*iCPs - 1) = dRI(iCPs,2);
                dRdxiMtxI(3,3*iCPs) = dRI(iCPs,2);

                % dR/deta
                dRdetaMtxI(1,3*iCPs - 2) = dRI(iCPs,3);
                dRdetaMtxI(2,3*iCPs - 1) = dRI(iCPs,3);
                dRdetaMtxI(3,3*iCPs) = dRI(iCPs,3);
            end
            
            % For patch J :
            % _____________
            
            for iCPs = 1:noCPsElJ
                % R
                RMtxJ(1,3*iCPs - 2) = dRJ(iCPs,1);
                RMtxJ(2,3*iCPs - 1) = dRJ(iCPs,1);
                RMtxJ(3,3*iCPs) = dRJ(iCPs,1);

                % dR/dxi
                dRdxiMtxJ(1,3*iCPs - 2) = dRJ(iCPs,2);
                dRdxiMtxJ(2,3*iCPs - 1) = dRJ(iCPs,2);
                dRdxiMtxJ(3,3*iCPs) = dRJ(iCPs,2);

                % dR/deta
                dRdetaMtxJ(1,3*iCPs - 2) = dRJ(iCPs,3);
                dRdetaMtxJ(2,3*iCPs - 1) = dRJ(iCPs,3);
                dRdetaMtxJ(3,3*iCPs) = dRJ(iCPs,3);
            end
            
            %% 4ii.6. Compute the covariant base vectors of the reference configuration
            
            % For patch I :
            % _____________
            
            [A1I,A2I] = computeBaseVectorsAndDerivativesForBSplineSurface...
                (xiSpanI,pI,etaSpanI,qI,CPI,0,dRI);
            A3TildeI = cross(A1I,A2I);
            A3I = A3TildeI/norm(cross(A1I,A2I));
            
            % For patch J :
            % _____________
            
            [A1J,A2J] = computeBaseVectorsAndDerivativesForBSplineSurface...
                (xiSpanJ,pJ,etaSpanJ,qJ,CPJ,0,dRJ);
            A3TildeJ = cross(A1J,A2J);
            A3J = A3TildeJ/norm(cross(A1J,A2J));
            
            
            %% 4ii.7. Compute the normal to the boundary vector and transform it into the contavariant basis
            
            % For patch I :
            % _____________
            
            [uGCI,~] = computeNormalAndTangentVectorsToBSplineBoundary...
                (xiI,XiI,etaI,EtaI,A1I,A2I,A3I,isOnXiI);
            uContravariantI = [A1I'
                               A2I']*uGCI;
                           
            % For patch J :
            % _____________
            
            [uGCJ,~] = computeNormalAndTangentVectorsToBSplineBoundary...
                (xiJ,XiJ,etaJ,EtaJ,A1J,A2J,A3J,isOnXiJ);
            uContravariantJ = [A1J'
                               A2J']*uGCJ;
                           
            %% 4ii.8. Compute the covariant base vectors of the current configuration
            
            % For patch I :
            % _____________
            
            [a1I,a2I] = computeBaseVectorsAndDerivativesForBSplineSurface...
                (xiSpanI,pI,etaSpanI,qI,CPdI,0,dRI);
            
            % For patch J :
            % _____________
            
            [a1J,a2J] = computeBaseVectorsAndDerivativesForBSplineSurface...
                (xiSpanJ,pJ,etaSpanJ,qJ,CPdJ,0,dRJ);
        
            %% 4ii.9. Compute the covariant metric coefficients of the current configuration
            
            % For patch I :
            % _____________
            
            aabCovariantI = [a1I a2I]'*[a1I a2I];
            
            % For patch J :
            % _____________
            
            aabCovariantJ = [a1J a2J]'*[a1J a2J];
            
            %% 4ii.10. Compute the contravariant bases
            
            % For patch I :
            % _____________
            
            AabCovariantI = [A1I A2I]'*[A1I A2I];
            AContravariantI = AabCovariantI\[A1I A2I]';
            AContravariantI = AContravariantI';
            
            % For patch J :
            % _____________
            
            AabCovariantJ = [A1J A2J]'*[A1J A2J];
            AContravariantJ = AabCovariantJ\[A1J A2J]';
            AContravariantJ = AContravariantJ';
            
            %% 4ii.11. Compute the local Cartesian bases
            
            % For patch I :
            % _____________
            
            eLCI = computeLocalCartesianBasis4BSplineSurface...
                ([A1I A2I],AContravariantI);
            
            % For patch J :
            % _____________
            
            eLCJ = computeLocalCartesianBasis4BSplineSurface...
                ([A1J A2J],AContravariantJ);
            
            %% 4ii.12. Compute the transformation matrix from the contravariant basis to the local Cartesian one
            
            % For patch I :
            % _____________
            
            TFromContraToLC4VoigtStrainI = ...
                computeTFromContra2LocalCartesian4VoigtStrainIGAKLShell...
                (eLCI,AContravariantI);
            
            % For patch J :
            % _____________
            
            TFromContraToLC4VoigtStrainJ = ...
                computeTFromContra2LocalCartesian4VoigtStrainIGAKLShell...
                (eLCJ,AContravariantJ);
            
            %% 4ii.13. Compute the transformation matrix from the local Cartesian bases to the covariant one
            
            % For patch I :
            % _____________
            
            TFromLCToCovI = computeTFromLocalCartesian2CovariantBasis4BSplineSurface...
                (eLCI,AContravariantI);
            
            % For patch I :
            % _____________
            
            TFromLCToCovJ = computeTFromLocalCartesian2CovariantBasis4BSplineSurface...
                (eLCJ,AContravariantJ);
            
            %% 4ii.14. Compute the Green-Lagrange strains in the contravariant bases
            
            % For patch I :
            % _____________
            
            EpsilonContraI = .5*[aabCovariantI(1,1) - AabCovariantI(1,1)
                                 aabCovariantI(2,2) - AabCovariantI(2,2)
                                 aabCovariantI(1,2) - AabCovariantI(1,2)];
                             
            % For patch I :
            % _____________
            
            EpsilonContraJ = .5*[aabCovariantJ(1,1) - AabCovariantJ(1,1)
                                 aabCovariantJ(2,2) - AabCovariantJ(2,2)
                                 aabCovariantJ(1,2) - AabCovariantJ(1,2)];
                             
            %% 4ii.15. Transform the Green-Lagrange strains in the local Cartesian bases
            
            % For patch I :
            % _____________
            
            EpsilonLCI = TFromContraToLC4VoigtStrainI*EpsilonContraI;
            
            % For patch I :
            % _____________
            
            EpsilonLCJ = TFromContraToLC4VoigtStrainJ*EpsilonContraJ;
            
            %% 4ii.16. Compute the prestress values on the local Cartesian coordinate systems
            
            % For patch I :
            % _____________
            
            % Check if a user defined coordinate system for the prestresses 
            % is chosen
            isPrestressOverDefinedSystemI = false;
            if isfield(parametersI.prestress,'computeBaseVectors')
                if ~isfield(parametersI.prestress,'computeParametricCoordinates')
                    error('Function handle parameters.prestress.computeParametricCoordinates has to be defined when defining the prestress over a user-defined coordinate system');
                end
                isPrestressOverDefinedSystemI = true;
            end
            
            % Compute the convective coordinates of the surface
            if isPrestressOverDefinedSystemI || isa(parametersI.prestress.voigtVector,'function_handle')
                XI = computeCartesianCoordinatesOfAPointOnBSplineSurface...
                    (xiSpanI,pI,xiI,XiI,etaSpanI,qI,etaI,EtaI,CPI,dRI(:,1));
                thetaI = prestressI.computeParametricCoordinates(XI);
            end
            
            % Compute the transformation matrix from the user defined 
            % coordinate system to the local Cartesian coordinate system if 
            % a user defined coordinate system is chosen
            if isPrestressOverDefinedSystemI
                prestressBaseVctI = prestressI.computeBaseVectors(thetaI(1,1),thetaI(2,1));
                T2LCI = computeT2LocalCartesianBasis(prestressBaseVctI,eLCI);
            else
                T2LCI = [1 0 0
                         0 1 0
                         0 0 1];
            end
            
            % Compute the prestress values
            if isa(prestressI.voigtVector,'function_handle')
                pTildeI = prestressI.voigtVector(thetaI);
            else
                pTildeI = prestressI.voigtVector;
            end
            
            % Transform the vector to the local Cartesian space if defined 
            % over a user defined coordinate system
            pTildeI = T2LCI*pTildeI;
            
            % For patch J :
            % _____________
            
            % Check if a user defined coordinate system for the prestresses 
            % is chosen
            isPrestressOverDefinedSystemJ = false;
            if isfield(parametersJ.prestress,'computeBaseVectors')
                if ~isfield(parametersJ.prestress,'computeParametricCoordinates')
                    error('Function handle parameters.prestress.computeParametricCoordinates has to be defined when defining the prestress over a user-defined coordinate system');
                end
                isPrestressOverDefinedSystemJ = true;
            end
            
            % Compute the convective coordinates of the surface
            if isPrestressOverDefinedSystemJ || isa(parametersJ.prestress.voigtVector,'function_handle')
                XJ = computeCartesianCoordinatesOfAPointOnBSplineSurface...
                    (xiSpanJ,pJ,xiJ,XiJ,etaSpanJ,qJ,etaJ,EtaJ,CPJ,dRJ(:,1));
                thetaJ = prestressJ.computeParametricCoordinates(XJ);
            end
            
            % Compute the transformation matrix from the user defined 
            % coordinate system to the local Cartesian coordinate system if 
            % a user defined coordinate system is chosen
            if isPrestressOverDefinedSystemJ
                prestressBaseVctJ = prestressJ.computeBaseVectors(thetaJ(1,1),thetaJ(2,1));
                T2LCJ = computeT2LocalCartesianBasis(prestressBaseVctJ,eLCJ);
            else
                T2LCJ = [1 0 0
                         0 1 0
                         0 0 1];
            end
            
            % Compute the prestress values
            if isa(prestressJ.voigtVector,'function_handle')
                pTildeJ = prestressJ.voigtVector(thetaJ);
            else
                pTildeJ = prestressJ.voigtVector;
            end
            
            % Transform the vector to the local Cartesian space if defined 
            % over a user defined coordinate system
            pTildeJ = T2LCJ*pTildeJ;
            
            %% 4ii.17. Compute the 2nd Piola-kirchhoff stresses in the local Cartesian systems
            
            % For patch I :
            % _____________
            
            NLCI = thicknessI*pTildeI + materialMtxVoigtI*EpsilonLCI;
            
            % For patch J :
            % _____________
            
            NLCJ = thicknessJ*pTildeJ +  materialMtxVoigtJ*EpsilonLCJ;
            
            %% 4ii.18. Transform the 2nd Piola-kirchhoff stresses in the covariant systems
            
            % For patch I :
            % _____________
            
            NCovariantI = TFromLCToCovI*NLCI;
            
            % For patch J :
            % _____________
            
            NCovariantJ = TFromLCToCovJ*NLCJ;
            
            %% 4ii.19. Compute the stress components
            
            % For patch I :
            % _____________
            
            PalphaBetaI = [NCovariantI(1,1) NCovariantI(3,1)
                           NCovariantI(3,1) NCovariantI(2,1)];
                       
            % For patch J :
            % _____________
            
            PalphaBetaJ = [NCovariantJ(1,1) NCovariantJ(3,1)
                           NCovariantJ(3,1) NCovariantJ(2,1)];
                       
            %% 4ii.20. Compute the traction vectors
            
            % For patch I :
            % _____________
            
            tractionVctI = [a1I a2I]*PalphaBetaI*uContravariantI;
            
            % For patch J :
            % _____________
            
            tractionVctJ = [a1J a2J]*PalphaBetaJ*uContravariantJ;
            
            %% 4ii.21. Compute the first variation of the Green-Lagrange strains in the contravariant bases with respect to the DOFs
            
            % For patch I :
            % _____________
            
            dEpsilonContraI = [a1I(:,1)'*dRdxiMtxI
                               a2I(:,1)'*dRdetaMtxI
                               .5*(a2I(:,1)'*dRdxiMtxI + a1I(:,1)'*dRdetaMtxI)];
                           
            % For patch J :
            % _____________
            
            dEpsilonContraJ = [a1J(:,1)'*dRdxiMtxJ
                               a2J(:,1)'*dRdetaMtxJ
                               .5*(a2J(:,1)'*dRdxiMtxJ + a1J(:,1)'*dRdetaMtxJ)];
                           
            %% 4ii.22. Transform the first variation of the Green-Lagrange strains at the local Cartesian bases
            
            % For patch I :
            % _____________
            
            dEpsilonCartesianI = TFromContraToLC4VoigtStrainI*dEpsilonContraI;
            
            % For patch J :
            % _____________
            
            dEpsilonCartesianJ = TFromContraToLC4VoigtStrainJ*dEpsilonContraJ;
            
            %% 4ii.23. Compute the first variations of the 2nd Piola-Kichhoff stresses in the local Cartesian bases
            
            % For patch I :
            % _____________
            
            dNCartesianI = materialMtxVoigtI*dEpsilonCartesianI;
            
            % For patch J :
            % _____________
            
            dNCartesianJ = materialMtxVoigtJ*dEpsilonCartesianJ;
            
            %% 4ii.24. Transform the first variations of the 2nd Piola-Kichhoff stresses at the covariant bases
            
            % For patch I :
            % _____________
            
            dNCovariantI = TFromLCToCovI*dNCartesianI;
            
            % For patch J :
            % _____________
            
            dNCovariantJ = TFromLCToCovJ*dNCartesianJ;
            
            %% 4ii.25. Compute the first variation of the traction vectors
            
            % For patch I :
            % _____________
            
            dtractionVctI = uContravariantI(1,1)*a1I*dNCovariantI(1,:) + ...
                uContravariantI(2,1)*a2I*dNCovariantI(2,:) + ...
                (uContravariantI(2,1)*a1I + uContravariantI(1,1)*a2I)*dNCovariantI(3,:) + ...
                (uContravariantI(1,1)*NCovariantI(1,1) + uContravariantI(2,1)*NCovariantI(3,1))*dRdxiMtxI + ...
                (uContravariantI(2,1)*NCovariantI(2,1) + uContravariantI(1,1)*NCovariantI(3,1))*dRdetaMtxI;
            
            % For patch J :
            % _____________
            
            dtractionVctJ = uContravariantJ(1,1)*a1J*dNCovariantJ(1,:) + ...
                uContravariantJ(2,1)*a2J*dNCovariantJ(2,:) + ...
                (uContravariantJ(2,1)*a1J + uContravariantJ(1,1)*a2J)*dNCovariantJ(3,:) + ...
                (uContravariantJ(1,1)*NCovariantJ(1,1) + uContravariantJ(2,1)*NCovariantJ(3,1))*dRdxiMtxJ + ...
                (uContravariantJ(2,1)*NCovariantJ(2,1) + uContravariantJ(1,1)*NCovariantJ(3,1))*dRdetaMtxJ;
            
            %% 4ii.26. Compute the necessary products needed for the second variations of the traction vectors
            
            % For patch I :
            % _____________
            
            PiMtxI = [a1I*uContravariantI(1,1) a2I*uContravariantI(2,1) a2I*uContravariantI(1,1) + a1I*uContravariantI(2,1)]*...
                TFromLCToCovI*materialMtxVoigtI*TFromContraToLC4VoigtStrainI;
            productsII = dispVctI'*PiMtxI;

            % For patch J :
            % _____________
            
            PiMtxJ = [a1J*uContravariantJ(1,1) a2J*uContravariantJ(2,1) a2J*uContravariantJ(1,1) + a1J*uContravariantJ(2,1)]*...
                TFromLCToCovJ*materialMtxVoigtJ*TFromContraToLC4VoigtStrainJ;
            productsJJ = dispVctJ'*PiMtxJ;
            
            % For coupling between patch I and J :
            % ____________________________________
            
            productsIJ = dispVctI'*PiMtxJ;
            productsJI = dispVctJ'*PiMtxI;
            
            %% 4ii.27. Compute the determinant of the Jacobian of the transformation from the physical to the parent domain on the GP
            if isOnXiI
                detJxxi = norm(A1I(:,1));
            else
                detJxxi = norm(A2I(:,1));
            end
            if detJxxi < tolDet
                continue;
            end
            
            %% 4ii.28. Compute the element length at the GP
            elementLengthOnGP = detJxxi*detJxizeta*GW(j);
            
            %% 4ii.29. Compute the element tangent stiffness matrix contributions on the GP and add them to the global matrices
            
            % For patch I :
            % _____________
            
            KNitscheI(EFTI,EFTI) = KNitscheI(EFTI,EFTI) - ...
                propCoupling.gammaTilde*(dtractionVctI'*RMtxI + RMtxI'*dtractionVctI + ...
                productsII(1,1)*(dRdxiMtxI'*dRdxiMtxI) + ...
                productsII(1,2)*(dRdetaMtxI'*dRdetaMtxI) + ...
                productsII(1,3)*.5*(dRdxiMtxI'*dRdetaMtxI + dRdetaMtxI'*dRdxiMtxI) + ...% [a1 a2]*ddNalphaBeta*uContravariant
                uContravariantI(1,1)*dRdxiMtxI'*dispVctI*dNCovariantI(1,:) + ...
                uContravariantI(2,1)*dRdetaMtxI'*dispVctI*dNCovariantI(2,:) + ...
                (uContravariantI(2,1)*dRdxiMtxI' + uContravariantI(1,1)*dRdetaMtxI')*dispVctI*dNCovariantI(3,:) + ... % daCovariant_s'*squeeze(dNalphaBeta(rDOFs,:,:))*uContravariant
                uContravariantI(1,1)*dNCovariantI(1,:)'*dispVctI'*dRdxiMtxI + ...
                uContravariantI(2,1)*dNCovariantI(2,:)'*dispVctI'*dRdetaMtxI + ...
                dNCovariantI(3,:)'*dispVctI'*(uContravariantI(2,1)*dRdxiMtxI + uContravariantI(1,1)*dRdetaMtxI) - ... % daCovariant_r'*squeeze(dNalphaBeta(sDOFs,:,:))*uContravariant
                productsJI(1,1)*(dRdxiMtxI'*dRdxiMtxI) - ...
                productsJI(1,2)*(dRdetaMtxI'*dRdetaMtxI) - ...
                productsJI(1,3)*.5*(dRdxiMtxI'*dRdetaMtxI + dRdetaMtxI'*dRdxiMtxI) - ... % [a1 a2]*ddNalphaBeta*uContravariant
                uContravariantI(1,1)*dRdxiMtxI'*dispVctJ*dNCovariantI(1,:) - ...
                uContravariantI(2,1)*dRdetaMtxI'*dispVctJ*dNCovariantI(2,:) - ...
                (uContravariantI(2,1)*dRdxiMtxI' + uContravariantI(1,1)*dRdetaMtxI')*dispVctJ*dNCovariantI(3,:) - ... % daCovariant_s'*squeeze(dNalphaBeta(rDOFs,:,:))*uContravariant
                uContravariantI(1,1)*dNCovariantI(1,:)'*dispVctJ'*dRdxiMtxI - ...
                uContravariantI(2,1)*dNCovariantI(2,:)'*dispVctJ'*dRdetaMtxI - ...
                dNCovariantI(3,:)'*dispVctJ'*(uContravariantI(2,1)*dRdxiMtxI + uContravariantI(1,1)*dRdetaMtxI))*elementLengthOnGP; ... % daCovariant_r'*squeeze(dNalphaBeta(sDOFs,:,:))*uContravariant
            
            % For patch J :
            % _____________
            
            KNitscheJ(EFTJ,EFTJ) = KNitscheJ(EFTJ,EFTJ) - ...
                propCoupling.gammaTilde*(dtractionVctJ'*RMtxJ + RMtxJ'*dtractionVctJ - ...
                productsIJ(1,1)*(dRdxiMtxJ'*dRdxiMtxJ) - ...
                productsIJ(1,2)*(dRdetaMtxJ'*dRdetaMtxJ) - ...
                productsIJ(1,3)*.5*(dRdxiMtxJ'*dRdetaMtxJ + dRdetaMtxJ'*dRdxiMtxJ) - ... % [a1 a2]*ddNalphaBeta*uContravariant
                uContravariantJ(1,1)*dRdxiMtxJ'*dispVctI*dNCovariantJ(1,:) - ...
                uContravariantJ(2,1)*dRdetaMtxJ'*dispVctI*dNCovariantJ(2,:) - ...
                (uContravariantJ(2,1)*dRdxiMtxJ' + uContravariantJ(1,1)*dRdetaMtxJ')*dispVctI*dNCovariantJ(3,:) - ... % daCovariant_s'*squeeze(dNalphaBeta(rDOFs,:,:))*uContravariant
                uContravariantJ(1,1)*dNCovariantJ(1,:)'*dispVctI'*dRdxiMtxJ - ...
                uContravariantJ(2,1)*dNCovariantJ(2,:)'*dispVctI'*dRdetaMtxJ - ...
                dNCovariantJ(3,:)'*dispVctI'*(uContravariantJ(2,1)*dRdxiMtxJ + uContravariantJ(1,1)*dRdetaMtxJ) + ...
                productsJJ(1,1)*(dRdxiMtxJ'*dRdxiMtxJ) + ...
                productsJJ(1,2)*(dRdetaMtxJ'*dRdetaMtxJ) + ...
                productsJJ(1,3)*.5*(dRdxiMtxJ'*dRdetaMtxJ + dRdetaMtxJ'*dRdxiMtxJ) + ... % [a1 a2]*ddNalphaBeta*uContravariant
                uContravariantJ(1,1)*dRdxiMtxJ'*dispVctJ*dNCovariantJ(1,:) + ...
                uContravariantJ(2,1)*dRdetaMtxJ'*dispVctJ*dNCovariantJ(2,:) + ...
                (uContravariantJ(2,1)*dRdxiMtxJ' + uContravariantJ(1,1)*dRdetaMtxJ')*dispVctJ*dNCovariantJ(3,:) + ... % daCovariant_s'*squeeze(dNalphaBeta(rDOFs,:,:))*uContravariant
                uContravariantJ(1,1)*dNCovariantJ(1,:)'*dispVctJ'*dRdxiMtxJ + ...
                uContravariantJ(2,1)*dNCovariantJ(2,:)'*dispVctJ'*dRdetaMtxJ + ...
                dNCovariantJ(3,:)'*dispVctJ'*(uContravariantJ(2,1)*dRdxiMtxJ + uContravariantJ(1,1)*dRdetaMtxJ))*elementLengthOnGP; % daCovariant_r'*squeeze(dNalphaBeta(sDOFs,:,:))*uContravariant
            
            %% 4ii.30. Compute the element tangent coupling matrices at the GP
        
            % For patch I :
            % _____________
                
            CNitscheI(EFTI,EFTJ) = CNitscheI(EFTI,EFTJ) - ...
                propCoupling.gammaTilde*(- dtractionVctI'*RMtxJ - RMtxI'*dtractionVctJ)*elementLengthOnGP;
            
            % For patch J :
            % _____________
            
            % CNitscheJ = CNitscheI' due to the symmetry of the variational
            % problem
            
            %% 4ii.31. Compute the element residual vector contributions on the GP and add the to the global vectors
            
            % For patch I :
            % _____________
            
            % resVctNitscheI(EFTI,1) = resVctNitscheI(EFTI,1) - 0.5*(RMtxI'*(tractionVctI - tractionVctJ) + dtractionVctI'*(dispVctI - dispVctJ))
            resVctNitscheI(EFTI,1) =  resVctNitscheI(EFTI,1) - ...
                propCoupling.gammaTilde*(RMtxI'*tractionVctI - RMtxI'*tractionVctJ + ...
                dtractionVctI'*dispVctI - dtractionVctI'*dispVctJ)*elementLengthOnGP;
            
            % For patch J :
            % _____________
            
            % resVctNitscheJ(EFTJ,1) = resVctNitscheJ(EFTJ,1) - 0.5*(RMtxJ'*(tractionVctJ - tractionVctI) + dtractionVctJ'*(dispVctJ - dispVctI))
            resVctNitscheJ(EFTJ,1) = resVctNitscheJ(EFTJ,1) - ...
                propCoupling.gammaTilde*(RMtxJ'*tractionVctJ - RMtxJ'*tractionVctI - ...
                dtractionVctJ'*dispVctI + dtractionVctJ'*dispVctJ)*elementLengthOnGP;
            
            %% 4ii.32. Compute the mass matrices necessary for the stabilization of the variational problem at the tangent stiffness matrices level
            if propCoupling.estimationStabilPrm == true
                
                % For patch I :
                % _____________
            
                massMtxKI(EFTI,EFTI) = massMtxKI(EFTI,EFTI) + ...
                    (RMtxI'*RMtxI)*elementLengthOnGP;
                
                % For patch J :
                % _____________
                
                massMtxKJ(EFTJ,EFTJ) = massMtxKJ(EFTJ,EFTJ) + ...
                    (RMtxJ'*RMtxJ)*elementLengthOnGP;
            end
            
            %% 4ii.33. Compute the mass matrices necessary for the stabilization of the variational problem at the tangent coupling matrices level
            
            % For patch I :
            % _____________
                
            if propCoupling.estimationStabilPrm == true
                massMtxCI(EFTI,EFTJ) = massMtxCI(EFTI,EFTJ) - ...
                    (RMtxI'*RMtxJ)*elementLengthOnGP;
            end
            
            % For patch J :
            % _____________
            
            % massMtxCJ = massMtxCI' due to the symmetry of the variational
            % problem
            
            %% 4ii.34. Compute the second variation of the traction vectors needed for the interface eigenvalue problem which needs to be solved for the estimation of the stabilization parameters
%             if noNonlinearIteration == 1 && propCoupling.estimationStabilPrm == true
%                 
%                 % For patch I :
%                 % _____________
%      
%                 % ddtractionVctI = dtractionVctI + ...
%                 %                 [dHatElI'*squeeze(ddtractionVctI(1,:,:))' 
%                 %                  dHatElI'*squeeze(ddtractionVctI(2,:,:))' 
%                 %                  dHatElI'*squeeze(ddtractionVctI(3,:,:))']
%                 ddtractionVctI = dtractionVctI + ...
%                     [dHatElI'*(PiMtxI(1,1)*(dRdxiMtxI'*dRdxiMtxI) + PiMtxI(1,2)*(dRdetaMtxI'*dRdetaMtxI) + PiMtxI(1,3)*.5*(dRdxiMtxI'*dRdetaMtxI + dRdetaMtxI'*dRdxiMtxI))
%                      dHatElI'*(PiMtxI(2,1)*(dRdxiMtxI'*dRdxiMtxI) + PiMtxI(2,2)*(dRdetaMtxI'*dRdetaMtxI) + PiMtxI(2,3)*.5*(dRdxiMtxI'*dRdetaMtxI + dRdetaMtxI'*dRdxiMtxI))
%                      dHatElI'*(PiMtxI(3,1)*(dRdxiMtxI'*dRdxiMtxI) + PiMtxI(3,2)*(dRdetaMtxI'*dRdetaMtxI) + PiMtxI(3,3)*.5*(dRdxiMtxI'*dRdetaMtxI + dRdetaMtxI'*dRdxiMtxI))] + ...
%                     [dHatElI'*(uContravariantI(1,1)*dRdxiMtxI(1,:)'*dNCovariantI(1,:) + uContravariantI(2,1)*dRdetaMtxI(1,:)'*dNCovariantI(2,:) + (uContravariantI(2,1)*dRdxiMtxI(1,:)' + uContravariantI(1,1)*dRdetaMtxI(1,:)')*dNCovariantI(3,:))
%                      dHatElI'*(uContravariantI(1,1)*dRdxiMtxI(2,:)'*dNCovariantI(1,:) + uContravariantI(2,1)*dRdetaMtxI(2,:)'*dNCovariantI(2,:) + (uContravariantI(2,1)*dRdxiMtxI(2,:)' + uContravariantI(1,1)*dRdetaMtxI(2,:)')*dNCovariantI(3,:))
%                      dHatElI'*(uContravariantI(1,1)*dRdxiMtxI(3,:)'*dNCovariantI(1,:) + uContravariantI(2,1)*dRdetaMtxI(3,:)'*dNCovariantI(2,:) + (uContravariantI(2,1)*dRdxiMtxI(3,:)' + uContravariantI(1,1)*dRdetaMtxI(3,:)')*dNCovariantI(3,:))] + ...
%                      [dHatElI'*(uContravariantI(1,1)*dNCovariantI(1,:)'*dRdxiMtxI(1,:) + uContravariantI(2,1)*dNCovariantI(2,:)'*dRdetaMtxI(1,:) + dNCovariantI(3,:)'*(uContravariantI(2,1)*dRdxiMtxI(1,:) + uContravariantI(1,1)*dRdetaMtxI(1,:)))
%                       dHatElI'*(uContravariantI(1,1)*dNCovariantI(1,:)'*dRdxiMtxI(2,:) + uContravariantI(2,1)*dNCovariantI(2,:)'*dRdetaMtxI(2,:) + dNCovariantI(3,:)'*(uContravariantI(2,1)*dRdxiMtxI(2,:) + uContravariantI(1,1)*dRdetaMtxI(2,:)))
%                       dHatElI'*(uContravariantI(1,1)*dNCovariantI(1,:)'*dRdxiMtxI(3,:) + uContravariantI(2,1)*dNCovariantI(2,:)'*dRdetaMtxI(3,:) + dNCovariantI(3,:)'*(uContravariantI(2,1)*dRdxiMtxI(3,:) + uContravariantI(1,1)*dRdetaMtxI(3,:)))];
%                     
%                 % For patch J :
%                 % _____________
% 
%                 % ddtractionVctJ = dtractionVctJ + ...
%                 %                 [dHatElJ'*squeeze(ddtractionVctJ(1,:,:))' 
%                 %                  dHatElJ'*squeeze(ddtractionVctJ(2,:,:))' 
%                 %                  dHatElJ'*squeeze(ddtractionVctJ(3,:,:))']
%                 ddtractionVctJ = dtractionVctJ + ...
%                     [dHatElJ'*(PiMtxJ(1,1)*(dRdxiMtxJ'*dRdxiMtxJ) + PiMtxJ(1,2)*(dRdetaMtxJ'*dRdetaMtxJ) + PiMtxJ(1,3)*.5*(dRdxiMtxJ'*dRdetaMtxJ + dRdetaMtxJ'*dRdxiMtxJ))
%                      dHatElJ'*(PiMtxJ(2,1)*(dRdxiMtxJ'*dRdxiMtxJ) + PiMtxJ(2,2)*(dRdetaMtxJ'*dRdetaMtxJ) + PiMtxJ(2,3)*.5*(dRdxiMtxJ'*dRdetaMtxJ + dRdetaMtxJ'*dRdxiMtxJ))
%                      dHatElJ'*(PiMtxJ(3,1)*(dRdxiMtxJ'*dRdxiMtxJ) + PiMtxJ(3,2)*(dRdetaMtxJ'*dRdetaMtxJ) + PiMtxJ(3,3)*.5*(dRdxiMtxJ'*dRdetaMtxJ + dRdetaMtxJ'*dRdxiMtxJ))] + ...
%                     [dHatElJ'*(uContravariantJ(1,1)*dRdxiMtxJ(1,:)'*dNCovariantJ(1,:) + uContravariantJ(2,1)*dRdetaMtxJ(1,:)'*dNCovariantJ(2,:) + (uContravariantJ(2,1)*dRdxiMtxJ(1,:)' + uContravariantJ(1,1)*dRdetaMtxJ(1,:)')*dNCovariantJ(3,:))
%                      dHatElJ'*(uContravariantJ(1,1)*dRdxiMtxJ(2,:)'*dNCovariantJ(1,:) + uContravariantJ(2,1)*dRdetaMtxJ(2,:)'*dNCovariantJ(2,:) + (uContravariantJ(2,1)*dRdxiMtxJ(2,:)' + uContravariantJ(1,1)*dRdetaMtxJ(2,:)')*dNCovariantJ(3,:))
%                      dHatElJ'*(uContravariantJ(1,1)*dRdxiMtxJ(3,:)'*dNCovariantJ(1,:) + uContravariantJ(2,1)*dRdetaMtxJ(3,:)'*dNCovariantJ(2,:) + (uContravariantJ(2,1)*dRdxiMtxJ(3,:)' + uContravariantJ(1,1)*dRdetaMtxJ(3,:)')*dNCovariantJ(3,:))] + ...
%                      [dHatElJ'*(uContravariantJ(1,1)*dNCovariantJ(1,:)'*dRdxiMtxJ(1,:) + uContravariantJ(2,1)*dNCovariantJ(2,:)'*dRdetaMtxJ(1,:) + dNCovariantJ(3,:)'*(uContravariantJ(2,1)*dRdxiMtxJ(1,:) + uContravariantJ(1,1)*dRdetaMtxJ(1,:)))
%                       dHatElJ'*(uContravariantJ(1,1)*dNCovariantJ(1,:)'*dRdxiMtxJ(2,:) + uContravariantJ(2,1)*dNCovariantJ(2,:)'*dRdetaMtxJ(2,:) + dNCovariantJ(3,:)'*(uContravariantJ(2,1)*dRdxiMtxJ(2,:) + uContravariantJ(1,1)*dRdetaMtxJ(2,:)))
%                       dHatElJ'*(uContravariantJ(1,1)*dNCovariantJ(1,:)'*dRdxiMtxJ(3,:) + uContravariantJ(2,1)*dNCovariantJ(2,:)'*dRdetaMtxJ(3,:) + dNCovariantJ(3,:)'*(uContravariantJ(2,1)*dRdxiMtxJ(3,:) + uContravariantJ(1,1)*dRdetaMtxJ(3,:)))];
%             end
            
            %% 4ii.35. Compute the Q matrix necessary for the eigenvalue problem related to the estimation of the stabilization parameter at the tangent stiffness matrices level
            if noNonlinearIteration == 1 && propCoupling.estimationStabilPrm == true && min([norm(A3TildeI) norm(A3TildeJ)]) >= tolSurfaceNormal
                % For patch I :
                % _____________
                
                % Compute necessary products
                productsQII = tractionVctI'*PiMtxI;
                productsQJI = tractionVctJ'*PiMtxI;
                
%                 QKMtxI(EFTI,EFTI) = QKMtxI(EFTI,EFTI) + ...
%                     (ddtractionVctI'*ddtractionVctI)*elementLengthOnGP;

                QKMtxI(EFTI,EFTI) = QKMtxI(EFTI,EFTI) + ...
                    (propCoupling.gammaTilde)^2*(2*(dtractionVctI'*dtractionVctI) + ... % dTraction'*dTraction
                    productsQII(1,1)*(dRdxiMtxI'*dRdxiMtxI) + ...
                    productsQII(1,2)*(dRdetaMtxI'*dRdetaMtxI) + ...
                    productsQII(1,3)*.5*(dRdxiMtxI'*dRdetaMtxI + dRdetaMtxI'*dRdxiMtxI) + ...% [a1 a2]*ddNalphaBeta*uContravariant
                    uContravariantI(1,1)*dRdxiMtxI'*tractionVctI*dNCovariantI(1,:) + ...
                    uContravariantI(2,1)*dRdetaMtxI'*tractionVctI*dNCovariantI(2,:) + ...
                    (uContravariantI(2,1)*dRdxiMtxI' + uContravariantI(1,1)*dRdetaMtxI')*tractionVctI*dNCovariantI(3,:) + ... % daCovariant_s'*squeeze(dNalphaBeta(rDOFs,:,:))*uContravariant
                    uContravariantI(1,1)*dNCovariantI(1,:)'*tractionVctI'*dRdxiMtxI + ...
                    uContravariantI(2,1)*dNCovariantI(2,:)'*tractionVctI'*dRdetaMtxI + ...
                    dNCovariantI(3,:)'*tractionVctI'*(uContravariantI(2,1)*dRdxiMtxI + uContravariantI(1,1)*dRdetaMtxI) - ... % daCovariant_r'*squeeze(dNalphaBeta(sDOFs,:,:))*uContravariant
                    productsQJI(1,1)*(dRdxiMtxI'*dRdxiMtxI) - ...
                    productsQJI(1,2)*(dRdetaMtxI'*dRdetaMtxI) - ...
                    productsQJI(1,3)*.5*(dRdxiMtxI'*dRdetaMtxI + dRdetaMtxI'*dRdxiMtxI) - ... % [a1 a2]*ddNalphaBeta*uContravariant
                    uContravariantI(1,1)*dRdxiMtxI'*tractionVctJ*dNCovariantI(1,:) - ...
                    uContravariantI(2,1)*dRdetaMtxI'*tractionVctJ*dNCovariantI(2,:) - ...
                    (uContravariantI(2,1)*dRdxiMtxI' + uContravariantI(1,1)*dRdetaMtxI')*tractionVctJ*dNCovariantI(3,:) - ... % daCovariant_s'*squeeze(dNalphaBeta(rDOFs,:,:))*uContravariant
                    uContravariantI(1,1)*dNCovariantI(1,:)'*tractionVctJ'*dRdxiMtxI - ...
                    uContravariantI(2,1)*dNCovariantI(2,:)'*tractionVctJ'*dRdetaMtxI - ...
                    dNCovariantI(3,:)'*tractionVctJ'*(uContravariantI(2,1)*dRdxiMtxI + uContravariantI(1,1)*dRdetaMtxI))*elementLengthOnGP; ... % daCovariant_r'*squeeze(dNalphaBeta(sDOFs,:,:))*uContravariant)*elementLengthOnGP;
                
                % For patch J :
                % _____________
                
                % Compute necessary products
                productsQIJ = tractionVctI'*PiMtxJ;
                productsQJJ = tractionVctJ'*PiMtxJ;
                
%                 QKMtxJ(EFTJ,EFTJ) = QKMtxJ(EFTJ,EFTJ) + ...
%                     (ddtractionVctJ'*ddtractionVctJ)*elementLengthOnGP;

                QKMtxJ(EFTJ,EFTJ) = QKMtxJ(EFTJ,EFTJ) + ...
                    (propCoupling.gammaTilde)^2*(2*(dtractionVctJ'*dtractionVctJ) - ...
                    productsQIJ(1,1)*(dRdxiMtxJ'*dRdxiMtxJ) - ...
                    productsQIJ(1,2)*(dRdetaMtxJ'*dRdetaMtxJ) - ...
                    productsQIJ(1,3)*.5*(dRdxiMtxJ'*dRdetaMtxJ + dRdetaMtxJ'*dRdxiMtxJ) - ... % [a1 a2]*ddNalphaBeta*uContravariant
                    uContravariantJ(1,1)*dRdxiMtxJ'*tractionVctI*dNCovariantJ(1,:) - ...
                    uContravariantJ(2,1)*dRdetaMtxJ'*tractionVctI*dNCovariantJ(2,:) - ...
                    (uContravariantJ(2,1)*dRdxiMtxJ' + uContravariantJ(1,1)*dRdetaMtxJ')*tractionVctI*dNCovariantJ(3,:) - ... % daCovariant_s'*squeeze(dNalphaBeta(rDOFs,:,:))*uContravariant
                    uContravariantJ(1,1)*dNCovariantJ(1,:)'*tractionVctI'*dRdxiMtxJ - ...
                    uContravariantJ(2,1)*dNCovariantJ(2,:)'*tractionVctI'*dRdetaMtxJ - ...
                    dNCovariantJ(3,:)'*tractionVctI'*(uContravariantJ(2,1)*dRdxiMtxJ + uContravariantJ(1,1)*dRdetaMtxJ) + ...
                    productsQJJ(1,1)*(dRdxiMtxJ'*dRdxiMtxJ) + ...
                    productsQJJ(1,2)*(dRdetaMtxJ'*dRdetaMtxJ) + ...
                    productsQJJ(1,3)*.5*(dRdxiMtxJ'*dRdetaMtxJ + dRdetaMtxJ'*dRdxiMtxJ) + ... % [a1 a2]*ddNalphaBeta*uContravariant
                    uContravariantJ(1,1)*dRdxiMtxJ'*tractionVctJ*dNCovariantJ(1,:) + ...
                    uContravariantJ(2,1)*dRdetaMtxJ'*tractionVctJ*dNCovariantJ(2,:) + ...
                    (uContravariantJ(2,1)*dRdxiMtxJ' + uContravariantJ(1,1)*dRdetaMtxJ')*tractionVctJ*dNCovariantJ(3,:) + ... % daCovariant_s'*squeeze(dNalphaBeta(rDOFs,:,:))*uContravariant
                    uContravariantJ(1,1)*dNCovariantJ(1,:)'*tractionVctJ'*dRdxiMtxJ + ...
                    uContravariantJ(2,1)*dNCovariantJ(2,:)'*tractionVctJ'*dRdetaMtxJ + ...
                    dNCovariantJ(3,:)'*tractionVctJ'*(uContravariantJ(2,1)*dRdxiMtxJ + uContravariantJ(1,1)*dRdetaMtxJ))*elementLengthOnGP; % daCovariant_r'*squeeze(dNalphaBeta(sDOFs,:,:))*uContravariant
            end
            
            %% 4ii.36. Compute the Q matrix necessary for the eigenvalue problem related to the estimation of the stabilization parameter at the tangent coupling matrices level
            if noNonlinearIteration == 1 && propCoupling.estimationStabilPrm == true
               
                % For patch I :
                % _____________
                                   
%                 QCMtxI(EFTI,EFTJ) = QCMtxI(EFTI,EFTJ) + ...
%                     (ddtractionVctI'*ddtractionVctJ)*elementLengthOnGP;

                QCMtxI(EFTI,EFTJ) = QCMtxI(EFTI,EFTJ) + ...
                    propCoupling.gammaTilde^2*(-2*dtractionVctI'*dtractionVctJ)*elementLengthOnGP;
                    
                % For patch J :
                % _____________
                
                % QCMtxJ = QCMtxI' due to the symmetry of the variational
                % problem
            end
        end
    end
end

%% 5. Find the DOFs on the current interface boundary if the interface eigenvalue problem for the estimation of the stabilization factor is to be solved
if noNonlinearIteration == 1 && propCoupling.estimationStabilPrm == true
    
    % For patch I :
    % _____________

    couplingDOFsI = [];
    for direction = 1:3
        couplingDOFsI = findDofs3D(couplingDOFsI,patchI.xicoup,patchI.etacoup,direction,CPI);
    end

    % For patch J :
    % _____________

    couplingDOFsJ = [];
    for direction = 1:3
        couplingDOFsJ = findDofs3D(couplingDOFsJ,patchJ.xicoup,patchJ.etacoup,direction,CPJ);
    end
end

%% 6. Estimate the stabilization parameter for the current patch coupling if the interface eigenvalue problem for the estimation of the stabilization factor is to be solved
if noNonlinearIteration == 1 && propCoupling.estimationStabilPrm == true
    eigenValues = eig([QKMtxI(couplingDOFsI,couplingDOFsI)  QCMtxI(couplingDOFsI,couplingDOFsJ)
                       QCMtxI(couplingDOFsI,couplingDOFsJ)' QKMtxJ(couplingDOFsJ,couplingDOFsJ)],...
                       [full(tanStiffMtxI(couplingDOFsI,couplingDOFsI))     zeros(length(couplingDOFsI),length(couplingDOFsJ))
                        zeros(length(couplingDOFsJ),length(couplingDOFsI))  full(tanStiffMtxJ(couplingDOFsJ,couplingDOFsJ))]);
    if ~isreal(eigenValues)
        eigenValues = real(eigenValues);
    end
	stabilPrm = 4*mI*max(eigenValues);
    if ~ischar(propTransientAnalysis)
        propCoupling.automaticStabilization(noConnection,noTimeStep) = stabilPrm;
    end
    if strcmp(outMsg,'outputEnabled')
        fprintf(strcat(tab,'Stabilization for patch connection %d estimated to %d \n'),noConnection,stabilPrm);
    end
    if noConnection == connections.No && strcmp(outMsg,'outputEnabled')
        fprintf('\n');
    end
elseif noNonlinearIteration ~= 1 && propCoupling.estimationStabilPrm == true
    if ~ischar(propTransientAnalysis)
        stabilPrm = propCoupling.automaticStabilization(noConnection,noTimeStep);
    end
end

%% 7. Stabilize accordingly the tangent and the coupling stiffness matrices if the interface eigenvalue problem for the estimation of the stabilization factor is to be solved
if propCoupling.estimationStabilPrm == true
    
    % For patch I :
    % _____________
    
    % Update the tangent stiffness matrix
    KNitscheI = KNitscheI + stabilPrm*massMtxKI;
    
    % Update the tangent coupling matrix
    CNitscheI = CNitscheI + stabilPrm*massMtxCI;
    
    % Update the residual vector
    resVctNitscheI = resVctNitscheI + stabilPrm*(massMtxKI*dHatI + massMtxCI*dHatJ);
    
    % For patch J :
    % _____________
    
    % Update the tangent stiffness matrix
    KNitscheJ = KNitscheJ + stabilPrm*massMtxKJ;
    
    % Update the tangent coupling matrix
    % CNitscheJ = CNitscheI' due to the symmetry of the variational problem
    
    % Update the residual vector
    resVctNitscheJ = resVctNitscheJ + stabilPrm*(massMtxCI'*dHatI + massMtxKJ*dHatJ);
end

end
