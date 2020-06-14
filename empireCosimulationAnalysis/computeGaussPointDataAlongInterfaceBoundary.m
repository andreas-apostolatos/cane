function [couplingCurveNumGP, trCurveMasterGPs, trCurveSlaveGPs, trCurveGPWeights, ...
    trCurveMasterGPTangents, trCurveSlaveGPTangents, trCurveGPJacobianProducts] = ...
    computeGaussPointDataAlongInterfaceBoundary ...
    (BSplinePatchI, BSplinePatchJ, xiExtensionI, etaExtensionI, ...
    xiExtensionJ, etaExtensionJ, isSameOrientation, propIntInterface, ...
    tab, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the Gauss Point data along a given interface boundary.
%
%                         Input :
%   BSplinePatchI,BSplinePatchJ : Structure arrays containing information 
%                                 on the master and slave B-Spline patches,
%                                       .p,q : The polynomial orders of the 
%                                              B-Spline surface in both 
%                                              parametric directions
%                                    .Xi,Eta : The knot vectors in both 
%                                              parametric directions
%                                        .CP : The set of control points 
%                                              and weights
%                                   .isNURBS : Flag on whether the basis is 
%                                              a NURBS or a B-Spline
%     xiExtensionI,xiExtensionJ : Extension along the xi-direction of the
%                                 master and slave inteface boundaries
%   etaExtensionI,etaExtensionJ : Extension along the eta-direction of the
%                                 master and slave inteface boundaries
%             isSameOrientation : Flag on whether the interface
%                                 parametrizations have the same 
%                                 orientation or not
%              propIntInterface : Integration properties,
%                                       .type : 'default', 'user'
%                                      .noGPs : Number of Gauss Points
%                           tab : Tabulation when printing message in the 
%                                 command window
%                        outMsg : Enables message outputting in the command 
%                                 window if chosen as 'outputEnabled'
%
%                        Output :
%            couplingCurveNumGP : Total number of Gauss points along the 
%                                 Dirichlet boundary
%              trCurveMasterGPs : Array containing the parametric locations 
%                                 of all GPs along the master interface
%                                 boundary
%               trCurveSlaveGPs : Array containing the parametric locations 
%                                 of all GPs along the slave interface
%                                 boundary
%              trCurveGPWeights : Array containing all the Gauss weights
%       trCurveMasterGPTangents : Array containing all the coordinates of 
%                                 the tangent vectors in the Cartesian 
%                                 space along the master interface boundary
%        trCurveSlaveGPTangents : Array containing all the coordinates of 
%                                 the tangent vectors in the Cartesian 
%                                 space along the master interface boundary
%     trCurveGPJacobianProducts : Array containing all the Jacobian 
%                                 products times the Gauss point weights at 
%                                 the Gauss point along the Dirichlet 
%                                 boundary
%
% Function layout :
%
% 0. Read input
%
% 1. Get the running and the fixed parameters on the patch boundary where the weak boundary conditions are applied
%
% 2. Issue the Gauss quadrature
%
% 3. Compute the total number of Gauss points and initialize output arrays
%
% 4. Loop over the elements on the parameter space where the weak boundary conditions are applied 
% ->
%    4i. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
%
%   4ii. Loop over the Gauss points
%   ->
%        4ii.1. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
%
%        4ii.2. Compute the NURBS basis functions
%
%        4ii.3. Compute the covariant base vectors of the reference configuration
%
%        4ii.4. Compute the determinant of the Jacobian of the transformation from the physical to the parent domain on the GP
%
%        4ii.5. Compute the element length at the GP
%
%        4ii.6. Update the Gauss Point counter
%   <-
% <-
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf([tab 'Creating Gauss Point data along interface boundary\n']);
end

%% 0. Read input

% For patch I :
% _____________

% Reassign the analysis arrays
pI = BSplinePatchI.p;
qI = BSplinePatchI.q;
XiI = BSplinePatchI.Xi;
EtaI = BSplinePatchI.Eta;
CPI = BSplinePatchI.CP;
isNURBSI = BSplinePatchI.isNURBS;
nxiI = length(CPI(:,1,1));
netaI = length(CPI(1,:,1));

% For patch J :
% _____________

pJ = BSplinePatchJ.p;
qJ = BSplinePatchJ.q;
XiJ = BSplinePatchJ.Xi;
EtaJ = BSplinePatchJ.Eta;
CPJ = BSplinePatchJ.CP;
isNURBSJ = BSplinePatchJ.isNURBS;
nxiJ = length(CPJ(:,1,1));
netaJ = length(CPJ(1,:,1));

% Define a basic tolerance value
tolBasic = 1e-2;

% Get a characteristic length for the interface
isOnXiI = false;
if etaExtensionI(1) == etaExtensionI(2)
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

% Tolerance for the determinant of the transformation from the parameter to
% the Cartesian space
tolDet = tolBasic*characteristicLength;

% Initialize counters
counterGPs = 1;

%% 1. Get the running and the fixed parameters on the patch interface and the coupling region

% For patch I :
% _____________

if etaExtensionI(1) == etaExtensionI(2)
    % Coupled region in xi-direction
    couplingRegionI = xiExtensionI;
    
    % Find the correct spans for the coupled region
    spanStartI = findKnotSpan(couplingRegionI(1),XiI,nxiI);
    spanEndI = findKnotSpan(couplingRegionI(2),XiI,nxiI) + 1;
    
    % Corresponding to the coupled region knot span
    couplingRegionOnKnotVectorI = XiI(spanStartI:spanEndI);
    
    % Fixed parameter on the parametric net
    etaI = etaExtensionI(1);
    
    % Find the span where xiEta it lies in
    etaSpanI = findKnotSpan(etaI,EtaI,netaI);
    
    % Flag on whether the coupling line is over xi
    isOnXiI = true;
else
    % Coupled region in eta-direction
    couplingRegionI = etaExtensionI;
    
    % Find the correct spans for the coupled region
    spanStartI = findKnotSpan(couplingRegionI(1),EtaI,netaI);   
    spanEndI = findKnotSpan(couplingRegionI(2),EtaI,netaI) + 1;   
    
    % Corresponding to the coupled region knot span
    couplingRegionOnKnotVectorI = EtaI(spanStartI:spanEndI);
    
    % Fixed parameter on the parametric net
    xiI = xiExtensionI(1);
    
    % Find the span where uv it lies in
    xiSpanI = findKnotSpan(xiI,XiI,nxiI);
    
    % Flag on whether the coupling line is over eta
    isOnXiI = false;
end

% For patch J :
% _____________

if etaExtensionJ(1) == etaExtensionJ(2)
	% Coupled region in xi-direction
    couplingRegionJ = xiExtensionJ;
    
    % Find the correct spans for the coupled region
    spanStartJ = findKnotSpan(couplingRegionJ(1),XiJ,nxiJ);   
    spanEndJ = findKnotSpan(couplingRegionJ(2),XiJ,nxiJ)+1; 
    
    % Corresponding to the coupled region knot span
    couplingRegionOnKnotVectorJ = XiJ(spanStartJ:spanEndJ);
    

    % Fixed parameter on the parametric net
    etaJ = etaExtensionJ(1);
    
    % Find the span where xiEta it lies in
    etaSpanJ = findKnotSpan(etaJ,EtaJ,netaJ);
    
    % Flag on whether the coupling line is over xi
    isOnXiJ = true;
else
    % Coupled region in eta-direction
    couplingRegionJ = etaExtensionJ;
    
    % Find the correct spans for the coupled region
    spanStartJ = findKnotSpan(couplingRegionJ(1),EtaJ,netaJ);   
    spanEndJ = findKnotSpan(couplingRegionJ(2),EtaJ,netaJ)+1;
    
    % Corresponding to the coupled region knot span
    couplingRegionOnKnotVectorJ = EtaJ(spanStartJ:spanEndJ);

    % Fixed parameter on the parametric net
    xiJ = xiExtensionJ(1);
    
    % Find the span where uv it lies in
    xiSpanJ = findKnotSpan(xiJ,XiJ,nxiJ);
    
    % Flag on whether the coupling line is over eta
    isOnXiJ = false;
end

% Merge the two knot vectors into one for integration purposes
couplingRegionOnKnotVector = mergesorted(couplingRegionOnKnotVectorI,couplingRegionOnKnotVectorJ);

% Delete double entries
couplingRegionOnKnotVector = unique(couplingRegionOnKnotVector);

%% 2. Issue the Gauss quadrature
if strcmp(propIntInterface.type,'default')
    if isOnXiI
        pDegreeI = pI;
    else
        pDegreeI = qI;
    end
    if isOnXiJ
        pDegreeJ = pJ;
    else
        pDegreeJ = qJ;
    end
    noGPs = ceil((pDegreeI + pDegreeJ + 1)/2);
elseif strcmp(propIntInterface.type,'user')
    noGPs = propIntInterface.noGPs;
end
[GP,GW] = getGaussPointsAndWeightsOverUnitDomain(noGPs);
if strcmp(outMsg,'outputEnabled')
    fprintf([tab '\t' 'Creating %d Gauss Points per knot span\n'],noGPs);
end

%% 3. Compute the total number of Gauss points and initialize output arrays
noGPsTotal = 0;
for i = 1:length(couplingRegionOnKnotVector)-1
    if couplingRegionOnKnotVector(i) ~= couplingRegionOnKnotVector(i+1)
        noGPsTotal = noGPsTotal + noGPs;
    end
end
if strcmp(outMsg,'outputEnabled')
    fprintf([tab '\t' 'Creating data for %d Gauss Points in total\n'],noGPsTotal);
end
couplingCurveNumGP = noGPsTotal;
trCurveMasterGPs = zeros(1,2*noGPsTotal);
trCurveSlaveGPs = zeros(1,2*noGPsTotal);
trCurveGPWeights = zeros(1,noGPsTotal);
trCurveMasterGPTangents = zeros(1,3*noGPsTotal);
trCurveSlaveGPTangents = zeros(1,3*noGPsTotal);
trCurveGPJacobianProducts = zeros(1,noGPsTotal);

%% 4. Loop over the elements on the parameter space where the weak boundary conditions are applied
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
            trCurveMasterGPs(1,2*counterGPs - 1) = xiI;
            trCurveMasterGPs(1,2*counterGPs) = etaI;
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
            trCurveSlaveGPs(1,2*counterGPs - 1) = xiJ;
            trCurveSlaveGPs(1,2*counterGPs) = etaJ;
            dRJ = computeIGABasisFunctionsAndDerivativesForSurface...
                (xiSpanJ,pJ,xiJ,XiJ,etaSpanJ,qJ,etaJ,EtaJ,CPJ,isNURBSJ,1);
                
            %% 4ii.3. Compute the covariant base vectors of the reference configuration
            
            % For patch I :
            % _____________
            
            [A1I,A2I] = computeBaseVectorsAndDerivativesForBSplineSurface...
                (xiSpanI,pI,etaSpanI,qI,CPI,0,dRI);
            if isOnXiI
                trCurveMasterGPTangents(1,3*counterGPs - 2) = A1I(1,1);
                trCurveMasterGPTangents(1,3*counterGPs - 1) = A1I(2,1);
                trCurveMasterGPTangents(1,3*counterGPs) = A1I(3,1);
            else
                trCurveMasterGPTangents(1,3*counterGPs - 2) = A2I(1,1);
                trCurveMasterGPTangents(1,3*counterGPs - 1) = A2I(2,1);
                trCurveMasterGPTangents(1,3*counterGPs) = A2I(3,1);
            end
            
            % For patch J :
            % _____________
            
            [A1J,A2J] = computeBaseVectorsAndDerivativesForBSplineSurface...
                (xiSpanJ,pJ,etaSpanJ,qJ,CPJ,0,dRJ);
            if isOnXiJ
                trCurveSlaveGPTangents(1,3*counterGPs - 2) = A1J(1,1);
                trCurveSlaveGPTangents(1,3*counterGPs - 1) = A1J(2,1);
                trCurveSlaveGPTangents(1,3*counterGPs) = A1J(3,1);
            else
                trCurveSlaveGPTangents(1,3*counterGPs - 2) = A2J(1,1);
                trCurveSlaveGPTangents(1,3*counterGPs - 1) = A2J(2,1);
                trCurveSlaveGPTangents(1,3*counterGPs) = A2J(3,1);
            end
            
            
            %% 4ii.4. Compute the determinant of the Jacobian of the transformation from the physical to the parent domain on the GP
            if isOnXiI
                detJxxi = norm(A1I(:,1));
            else
                detJxxi = norm(A2I(:,1));
            end
            if detJxxi < tolDet
                warning('Base vector with almost zero length was deteceted');
            end
            
            %% 4ii.5. Compute the element length at the GP
            elementLengthOnGP = detJxxi*detJxizeta*GW(j);
            trCurveGPWeights(1,counterGPs) = GW(j);
            trCurveGPJacobianProducts(1,counterGPs) = elementLengthOnGP;
            
            %% 4ii.6. Update the Gauss Point counter
            counterGPs = counterGPs + 1;
        end
    end
end

end
