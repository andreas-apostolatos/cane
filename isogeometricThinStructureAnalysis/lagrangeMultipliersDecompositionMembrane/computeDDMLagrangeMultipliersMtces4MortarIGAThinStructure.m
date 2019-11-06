%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _______________________________________________________               %
%   _______________________________________________________               %
%                                                                         %
%   Technische Universität München                                        %
%   Lehrstuhl für Statik, Prof. Dr.-Ing. Kai-Uwe Bletzinger               %
%   _______________________________________________________               %
%   _______________________________________________________               %
%                                                                         %
%                                                                         %
%   Authors                                                               %
%   _______________________________________________________________       %
%                                                                         %
%   Dipl.-Math. Andreas Apostolatos    (andreas.apostolatos@tum.de)       %
%   Dr.-Ing. Roland Wüchner            (wuechner@tum.de)                  %
%   Prof. Dr.-Ing. Kai-Uwe Bletzinger  (kub@tum.de)                       %
%   _______________________________________________________________       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [LambdaMaster,LambdaSlave] = computeDDMLagrangeMultipliersMtces4MortarIGAThinStructure...
    (patchMaster,patchSlave,haveSameOrientation,propCoupling)
%% Function documentation
%
% Returns the Lagrange Multipliers matrices corresponding to the 
% application of the mortar method for the multipatch coupling of
% thin-walled isogeometric structures.
%
%          Input :
% patchI, patchJ : 
%
%
%% Function main body

%% 0. Read input

% Master :
% --------

% Re-assign the patch properties
pMaster = patchMaster.p;
qMaster = patchMaster.q;
XiMaster = patchMaster.Xi;
EtaMaster = patchMaster.Eta;
CPMaster = patchMaster.CP;
isNURBSMaster = patchMaster.isNURBS;
xicoupMaster = patchMaster.xicoup;
etacoupMaster = patchMaster.etacoup;
noDOFsMaster = patchMaster.noDOFs;

% Number of knots in xi-,-eta directions
% mxiMaster = length(XiMaster);
% metaMaster = length(EtaMaster);

% Number of Control Points in xi-,eta- directions
nxiMaster = length(CPMaster(:,1,1));
netaMaster = length(CPMaster(1,:,1));

% Local number of control points that influence the element for the
% membrane surface
noCPsLocMaster = (pMaster + 1)*(qMaster + 1);

% Local number of degrees of freedom for the membrane surface
noDOFsLocMaster = 3*noCPsLocMaster;

% Get the patch numbering
DOFNumberingMaster = patchMaster.DOFNumbering;

% Initialize the element freedom table
EFTMaster = zeros(1,noDOFsLocMaster);

% Initialize auxiliary arrays
RMatrixMaster = zeros(3,noDOFsLocMaster);

% Slave : 
% -------

% Re-assign the patch properties
pSlave = patchSlave.p;
qSlave = patchSlave.q;
XiSlave = patchSlave.Xi;
EtaSlave = patchSlave.Eta;
CPSlave = patchSlave.CP;
isNURBSSlave = patchSlave.isNURBS;
xicoupSlave = patchSlave.xicoup;
etacoupSlave = patchSlave.etacoup;
noDOFsSlave = patchSlave.noDOFs;

% Number of knots in xi-,-eta directions
mxiSlave = length(XiSlave);
metaSlave = length(EtaSlave);

% Number of Control Points in xi-,eta- directions
nxiSlave = length(CPSlave(:,1,1));
netaSlave = length(CPSlave(1,:,1));

% Local number of control points that influence the element for the
% membrane surface
noCPsLocSlave = (pSlave + 1)*(qSlave + 1);

% Local number of degrees of freedom for the membrane surface
noDOFsLocSlave = 3*noCPsLocSlave;

% Get the patch numbering
DOFNumberingSlave = patchSlave.DOFNumbering;

% Initialize the element freedom table
EFTSlave = zeros(1,noDOFsLocSlave);

% Initialize auxiliary arrays
RMatrixSlave = zeros(3,noDOFsLocSlave);

% Define a singularity tolerance
tolArea = 1e-8;

% Initialize output arrays
LambdaMaster = zeros(noDOFsMaster,noDOFsSlave);
LambdaSlave = zeros(noDOFsSlave,noDOFsSlave);

%% 1. Get the running and the fixed parameters on the patch interface and the coupling region

% Master :
% --------

if etacoupMaster(1) == etacoupMaster(2)
    % Coupled region in xi-direction
    couplingRegionMaster = xicoupMaster;
    
    % Find the correct spans for the coupled region
    spanStartMaster = findKnotSpan(couplingRegionMaster(1),XiMaster,nxiMaster);
    spanEndMaster = findKnotSpan(couplingRegionMaster(2),XiMaster,nxiMaster) + 1;
    
    % Corresponding to the coupled region knot span
    couplingRegionOnKnotVectorMaster = XiMaster(spanStartMaster:spanEndMaster);
    
    % Fixed parameter on the parametric net
    etaMaster = etacoupMaster(1);
    
    % Find the span where xiEta it lies in
    etaSpanMaster = findKnotSpan(etaMaster,EtaMaster,netaMaster);
    
    % Flag on whether the coupling line is over xi
    isOnXiMaster = true;
else
    % Coupled region in eta-direction
    couplingRegionMaster = etacoupMaster;
    
    % Find the correct spans for the coupled region
    spanStartMaster = findKnotSpan(couplingRegionMaster(1),EtaMaster,netaMaster);   
    spanEndMaster = findKnotSpan(couplingRegionMaster(2),EtaMaster,netaMaster) + 1;
    
    % Corresponding to the coupled region knot span
    couplingRegionOnKnotVectorMaster = EtaMaster(spanStartMaster:spanEndMaster);
    
    % Fixed parameter on the parametric net
    xiMaster = xicoupMaster(1);
    
    % Find the span where xiEta it lies in
    xiSpanMaster = findKnotSpan(xiMaster,XiMaster,nxiMaster);
    
    % Flag on whether the coupling line is over eta
    isOnXiMaster = false;
end

% Slave :
% -------

if etacoupSlave(1) == etacoupSlave(2)
    % Coupled region in xi-direction
%     couplingRegionSlave = xicoupSlave;
    
    % Find the correct spans for the coupled region
%     spanStartSlave = findKnotSpan(couplingRegionSlave(1),XiSlave,nxiSlave);
%     spanEndSlave = findKnotSpan(couplingRegionSlave(2),XiSlave,nxiSlave) + 1;
    
    % Corresponding to the coupled region knot span
%     couplingRegionOnKnotVectorSlave = XiMaster(spanStartSlave:spanEndSlave);
    
    % Fixed parameter on the parametric net
    etaSlave = etacoupSlave(1);
    
    % Find the span where xiEta it lies in
    etaSpanSlave = findKnotSpan(etaSlave,EtaSlave,netaSlave);
    
    % Flag on whether the coupling line is over xi
    isOnXiSlave = true;
else
    % Coupled region in eta-direction
%     couplingRegionSlave = etacoupSlave;
    
    % Find the correct spans for the coupled region
%     spanStartSlave = findKnotSpan(couplingRegionSlave(1),EtaSlave,netaSlave);
%     spanEndSlave = findKnotSpan(couplingRegionSlave(2),EtaSlave,netaSlave) + 1;
    
    % Corresponding to the coupled region knot span
%     couplingRegionOnKnotVectorSlave = EtaSlave(spanStartSlave:spanEndSlave);
    
    % Fixed parameter on the parametric net
    xiSlave = xicoupSlave(1);
    
    % Find the span where xiEta it lies in
    xiSpanSlave = findKnotSpan(xiSlave,XiSlave,nxiSlave);
    
    % Flag on whether the coupling line is over eta
    isOnXiSlave = false;
end

%% 2. Compute the merged knot vector from the master and the slave patch

% Projected knot vector for the Lagrange multipliers field on the
% structural coupled field:
% If UL = [v1 v2] and US = [u1 u2] the two knot vectors, then the
% transformation rule is t(v) = u1 + (u1-u2)*(v-v1)/(v1-v2)

% Initialization of the projected knot vector
if isOnXiSlave
    knotVctSlaveProjected = XiSlave;
else
    knotVctSlaveProjected = EtaSlave;
end

% Assign the new entries of the projected Lagrange multipliers knot vector
for i = 1:length(XiSlave)
   knotVctSlaveProjected(i) = couplingRegionMaster(1) + (couplingRegionMaster(1)-couplingRegionMaster(2))*(XiSlave(i)-XiSlave(1))/(XiSlave(1)-XiSlave(end));
end

% Coupled parameter knot vector for the Lagrange multipliers field
if isOnXiSlave
    coupledRegionSlave = knotVctSlaveProjected(pSlave + 2:mxiSlave - pSlave - 1);
else
    coupledRegionSlave = knotVctSlaveProjected(qSlave + 2:metaSlave - qSlave - 1);
end

% Coupled parameter knot vector for the integration
couplingRegionOnKnotVectorTemp = mergesorted(couplingRegionOnKnotVectorMaster,coupledRegionSlave);
couplingRegionOnKnotVector = unique(couplingRegionOnKnotVectorTemp);

%% 3. Issue Gauss Point coordinates and weights
if strcmp(propCoupling.intC.type,'default')
    if isOnXi
        pDegree = p + 1;
    else
        pDegree = q + 1;
    end
    pLM = pLambda;
    if isLagrangeMultipliersRotationsEnabled
        pLM = max(pLM,pMu);
    end
    noGPs = ceil((pDegree + pLM + 1)/2);
elseif strcmp(propCoupling.intC.type,'user')
    noGPs = propCoupling.intC.noGPs;
end
[GP,GW] = getGaussPointsAndWeightsOverUnitDomain(noGPs);

%% 4. Loop over all the elements on the coupling interface
for i = 1:length(couplingRegionOnKnotVector) - 1
    %% 4i. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
    detJxizeta = (couplingRegionOnKnotVector(i + 1) - couplingRegionOnKnotVector(i))/2;
    
    %% 4ii. Loop over all Gauss points
    for j = 1:length(GP)
        %% 4ii.1. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
        xiEta = ((1-GP(j))*couplingRegionOnKnotVector(i)+(1+GP(j))*couplingRegionOnKnotVector(i+1))/2;
        
        %% 4ii.2. Compute the NURBS basis functions and their derivatives
        
        % Master :
        % --------
        
        if isOnXiMaster
            xiMaster = xiEta;
            xiSpanMaster = findKnotSpan(xiMaster,XiMaster,nxiMaster);
        else
            etaMaster = xiEta;
            etaSpanMaster = findKnotSpan(etaMaster,EtaMaster,netaMaster);
        end
        noDrvsBasis = 1;
        dRMaster = computeIGABasisFunctionsAndDerivativesForSurface...
            (xiSpanMaster,pMaster,xiMaster,XiMaster,etaSpanMaster,qMaster,etaMaster,EtaMaster,CPMaster,isNURBSMaster,noDrvsBasis);
        
        % Slave :
        % -------
        
        if ~haveSameOrientation
            if isOnXiSlave
                xiSlave = xiEta;
                xiSpanSlave = findKnotSpan(xiSlave,XiSlave,nxiSlave);
            else
                etaSlave = xiEta;
                etaSpanSlave = findKnotSpan(etaSlave,EtaSlave,netaSlave);
            end
        else
            if isOnXiSlave
                xiSlave = XiSlave(end) - xiEta;
                xiSpanSlave = findKnotSpan(xiSlave,XiSlave,nxiSlave);
            else
                etaSlave = EtaSlave(end) - xiEta;
                etaSpanSlave = findKnotSpan(etaSlave,EtaSlave,netaSlave);
            end
        end
        noDrvsBasis = 1;
        dRSlave = computeIGABasisFunctionsAndDerivativesForSurface...
            (xiSpanSlave,pSlave,xiSlave,XiSlave,etaSpanSlave,qSlave,etaSlave,EtaSlave,CPSlave,isNURBSSlave,noDrvsBasis);
        
        %% 4ii.3. Create the element freedom tables
        
        % Master :
        % --------

        % Initialize of the counter
        r = 1;

        % Relation global-local DoFs
        for cpj = etaSpanMaster - qMaster:etaSpanMaster
            for cpi = xiSpanMaster - pMaster:xiSpanMaster
                EFTMaster(r) = DOFNumberingMaster(cpi,cpj,1);
                EFTMaster(r + 1) = DOFNumberingMaster(cpi,cpj,2);
                EFTMaster(r + 2) = DOFNumberingMaster(cpi,cpj,3);

                % update counter
                r = r + 3;
            end
        end
        
        % Master :
        % --------

        % Initialize of the counter
        r = 1;

        % Relation global-local DoFs
        for cpj = etaSpanSlave - qSlave:etaSpanSlave
            for cpi = xiSpanSlave - pSlave:xiSpanSlave
                EFTSlave(r) = DOFNumberingSlave(cpi,cpj,1);
                EFTSlave(r + 1) = DOFNumberingSlave(cpi,cpj,2);
                EFTSlave(r + 2) = DOFNumberingSlave(cpi,cpj,3);

                % update counter
                r = r + 3;
            end
        end
        
        %% 4ii.4. Compute the covariant base vectors from the master side
        noDrvsBaseVct = 0;
        [dA1,dA2] = computeBaseVectorsAndDerivativesForBSplineSurface...
            (xiSpanMaster,pMaster,etaSpanMaster,qMaster,CPMaster,noDrvsBaseVct,dRMaster);
        if norm(cross(dA1(:,1),dA2(:,1))) < tolArea
            continue;
        end
        
        %% 4ii.5. Compute the basis functions matrices
        
        % Master :
        % --------

        k = 0;
        for c = 0:qMaster
            for b = 0:pMaster
                k = k + 1;
                RMatrixMaster(1,3*k - 2) = dRMaster(k,1);
                RMatrixMaster(2,3*k - 1) = dRMaster(k,1);
                RMatrixMaster(3,3*k) = dRMaster(k,1);
            end
        end
        
        % Slave :
        % --------

        k = 0;
        for c = 0:qSlave
            for b = 0:pSlave
                k = k + 1;
                RMatrixSlave(1,3*k - 2) = dRSlave(k,1);
                RMatrixSlave(2,3*k - 1) = dRSlave(k,1);
                RMatrixSlave(3,3*k) = dRSlave(k,1);
            end
        end
        
        %% 4ii.6. Compute the determinant of the Jacobian of the transformation from the physical to the parent domain on the GP
        if isOnXiMaster
            detJxxi = norm(dA1(:,1));
        else
            detJxxi = norm(dA2(:,1));
        end
        
        %% 4ii.7. Compute the element length at the GP
        elementLengthOnGP = detJxxi*detJxizeta*GW(j);
        
        %% 4ii.8. Compute the Lagrange Multipliers matrices and add them to the global matrices
        
        % Master :
        % --------
        
        LambdaMaster(EFTMaster,EFTSlave) = LambdaMaster(EFTMaster,EFTSlave) + ...
            RMatrixMaster'*RMatrixSlave*elementLengthOnGP;
        
        % Slave :
        % -------
        
        LambdaSlave(EFTSlave,EFTSlave) = LambdaSlave(EFTSlave,EFTSlave) - ...
            RMatrixSlave'*RMatrixSlave*elementLengthOnGP;
    end
end

end