function [relErrorL2GeoDomain,relErrorL2GeoInterface] = ...
    computeDomainAndInterfaceErrorInL2NormMembraneFormFiding...
    (BSplinePatches,connections,propReferenceSolution,propNewtonRapshon,...
    propError,propInt,outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the relative geometry error for each patch and the relative error
% for each interface between the patches with respect to the solution of a
% single patch problem considered as reference for the form finding of a
% multipatch membrane.
%
%                  Input :
%         BSplinePatches : Array of B-Spline surface patches
%            connections : Define the connection between the patches:
%                            .No : Number of connections
%                     .xiEtaCoup : [patchID1 patchID2 xi12 eta12 xi21 eta21
%                                   ...      ...    ...   ...  ...   ...]
%  propReferenceSolution : Properties of the reference solution. This can
%                          be defined either in terms of a patch which is
%                          highly refined or as an analytical function:
%                           .referenceBSplinePatch : A reference B-Spline
%                                                    patch or 'undefined'
%                          .computeRerenceSolution : Function handle for
%                                                    the computation of a
%                                                    reference solution
%                          solution is computed. If this is chosen to be 
%                          'undefined' no such error measure is computed
%      propNewtonRapshon : Properties for the Newton-Raphson algorithm 
%                          employed : 
%                            .eps : Residual tolerance
%                          .maxIt : Maximum number of iterations
%              propError : Properties on the computation of the error :
%                            .noSamplingPoints : Number of sampling points
%                                                for the initial guesses in
%                                                order to find the
%                                                corresponding Gauss point
%                                                on the reference B-Spline
%                                                patch
%                                    .tolClose : Tolerance up to which two
%                                                points are considered to 
%                                                be close, needed for the
%                                                Newton-Rapshon method for
%                                                the projection of the 
%                                                Gauss points
%                propInt : On the numerical integration over the domain :
%                            .type : 'default' or 'user'
%                         .noGPsXi : Number of Gauss points in xi-direction
%                        .noGPsEta : Number of Gauss points in eta-
%                                    direction
%                 outMsg : Enables outputting information onto the command
%                          window when chosen as 'outputEnabled'
%
%                 Output :
%    relErrorL2GeoDomain : Array containing the relative domain error for
%                          each patch
% relErrorL2GeoInterface : Array containing the interface error for each
%                          interface
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the patches in the multipatch geometry
% ->
%    1i. Get the patch parameters
%
%   1ii. Choose the Gauss Integration rule
%
%  1iii. Loop over the elements
%  ->
%        1iii.1. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
%
%        1iii.2. Loop over the Gauss points
%        ->
%                1iii.2i. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
%
%               1iii.2ii. Find the knot span
%
%              1iii.2iii. Compute the NURBS basis functions for the patch
%
%               1iii.2iv. Compute the physical location of the Gauss Point on patch
%
%                1iii.2v. Compute the parameter location of the Gauss Point on patch in patchReference
%
%               1iii.2vi. Find the knot span indices for the reference patch
%
%              1iii.2vii. Compute the NURBS basis functions for the reference patch
%
%             1iii.2viii. Compute the base vectors at the Gauss points for the patch
%
%               1iii.2ix. Compute the Cartesian coordinates of the Gauss points
%
%                1iii.2x. Compute the element area at the GP
%
%               1iii.2xi. Compute the L2-norms
%        <-
%  <-
% <-
%
% 2. Loop over all the interfaces in the multipatch geometry
% ->
%    2i. Get the patch IDs
%
%   2ii. Get the patch parameters
%
%  2iii. Check if the patches have the same or opposite orientation
%
%   2iv. Select an integration scheme
%
%    2v. Get the running and the fixed parameters on the patch interface and the coupling region
%
%   2vi. Compute the merged knot vector from both patches over the interface
%
%  2vii. Loop over the elements on the coupling surface
%  ->
%        2vii.1. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
%
%        2vii.2. Loop over the Gauss points
%        ->
%                2vii.2i. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
%
%               2vii.2ii. Compute the NURBS basis functions
%
%              2vii.2iii. Compute the edge length on the Gauss Point
%
%               2vii.2iv. Compute the physical location of the Gauss point
%
%                2vii.2v. Compute the parameter location of the Gauss Point on patch in patchReference
%
%               2iii.2vi. Find the knot span indices for the reference patch
%
%              2iii.2vii. Compute the NURBS basis functions for the reference patch
%
%             2vii.2viii. Compute the Cartesian coordinates of the Gauss points
%
%               2vii.2ix. Compute the element length at the GP
%
%                2vii.2x. Compute the L2-norms
%        <-
%  <-
% <-
%
% 3. Get the relative values
%
% 4. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('_____________________________________________________________________\n');
    fprintf('#####################################################################\n');
    fprintf('Computation of the constituents of the broken norm corresponding to \n');
    fprintf('the form-finding analysis of a multipatch membrane has been initiated \n');
    fprintf('\n');
    if ~isfield(propError,'noSamplingPoints')
        error('Define noSamplingPoints in the propError array');
    end
    fprintf('Number of sampling points which are chosen as initial guesses for the\n');
    fprintf('projection of the Gauss points is chosen as %d\n',propError.noSamplingPoints);
    if ~isfield(propError,'tolClose')
        error('Define tolClose in the propError array');
    end
    fprintf('Tolerance for point projection is chosen as %d\n',propError.tolClose);
    fprintf('_____________________________________________________________________\n');
    fprintf('\n');
    
    % start measuring computational time
    tic;
end

%% 0. Read input

% Number of B-Spline patches
noPatches = length(BSplinePatches);

% number of connections
noConnections = connections.No;

% Check input
isReferenceBSplinePatchProvided = false;
if isfield(propReferenceSolution,'referenceBSplinePatch')
    if ~strcmp(propReferenceSolution.referenceBSplinePatch,'undefined')
        isReferenceBSplinePatchProvided = true;
        referenceBSplinePatch = propReferenceSolution.referenceBSplinePatch;
    end
end
isReferenceSolutionProvided = false;
if isfield(propReferenceSolution,'computeRerenceSolution')
    if isa(propReferenceSolution.computeRerenceSolution,'function_handle')
        isReferenceSolutionProvided = true;
    end
end
if isReferenceBSplinePatchProvided && isReferenceSolutionProvided
    error('Either a reference patch or an analytical solution should be provided, not both simultaneously');
elseif ~isReferenceBSplinePatchProvided && ~isReferenceSolutionProvided
    error('Neither a reference patch nor an analytical solution is be provided');
end

% For the reference solution :
% ____________________________

if isReferenceBSplinePatchProvided
    % Recover the data from the reference B-Spline patch
    pRef = referenceBSplinePatch.p;
    qRef = referenceBSplinePatch.q;
    XiRef = referenceBSplinePatch.Xi;
    EtaRef = referenceBSplinePatch.Eta;
    CPRef = referenceBSplinePatch.CP;
    isNURBSRef = referenceBSplinePatch.isNURBS;

    % Number of Control Points in xi-,eta- directions
    nxiRef = length(CPRef(:,1,1));
    netaRef = length(CPRef(1,:,1));
else
    if isReferenceSolutionProvided
        computeReferenceSolution = propReferenceSolution.computeReferenceSolution;
    else
        error('Define computation of reference solution in propReferenceSolution.computeRerenceSolution as function handle');
    end
end

% Define a tolerance
tol_null = 1e-12;

% On the Netwon-Raphson iteration parameters
xi0Init = .5;
eta0Init = .5;

% On the initial guess for the Newton-Raphson iterations
dxi0 = (XiRef(end) - XiRef(1))/ceil(sqrt(propError.noSamplingPoints));
deta0 = (EtaRef(end) - EtaRef(1))/ceil(sqrt(propError.noSamplingPoints));

% Initialize the reference values for the displacement vector, the force 
% and the moment tensor if those are not provided

% Initialize output values
diffGeoDomainL2 = zeros(noPatches,1);
geoRefDomainL2 = zeros(noPatches,1);
diffGeoInterfaceL2 = zeros(connections.No,1);
geoInterfaceL2 = zeros(connections.No,1);

%% 1. Loop over all the patches in the multipatch geometry
for iPatches = 1:noPatches
    %% 1i. Get the patch parameters
    p = BSplinePatches{iPatches}.p;
    q = BSplinePatches{iPatches}.q;
    Xi = BSplinePatches{iPatches}.Xi;
    Eta = BSplinePatches{iPatches}.Eta;
    CP = BSplinePatches{iPatches}.CP;
    isNURBS = BSplinePatches{iPatches}.isNURBS;
    mxi = length(Xi);
    meta = length(Eta);
    nxi = length(CP(:,1,1));
    neta = length(CP(1,:,1));
    
    %% 1ii. Choose the Gauss Integration rule
    if strcmp(propInt.type,'default')
        noGPsXi = p + 1;
        noGPsEta = q + 1;
    elseif strcmp(propInt.type,'user')
        noGPsXi = propInt.noGPXiError;
        noGPsEta = propInt.noGPEtaError;
    end
    [GPXi,GWXi] = getGaussPointsAndWeightsOverUnitDomain(noGPsXi);
    [GPEta,GWEta] = getGaussPointsAndWeightsOverUnitDomain(noGPsEta);

    %% 1iii. Loop over the elements
    for iEta = q+1:meta-q-1
        for iXi = p+1:mxi-p-1
            if Xi(iXi+1) ~= Xi(iXi) && Eta(iEta+1) ~= Eta(iEta)
                %% 1iii.1. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
                detJxizeta = (Xi(iXi+1)-Xi(iXi))*(Eta(iEta+1)-Eta(iEta))/4;

                %% 1iii.2. Loop over the Gauss points
                for iGPEta = 1:length(GWEta)
                    for iGPXi = 1:length(GWXi)
                        %% 1iii.2i. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
                        xi = ( Xi(iXi+1)+Xi(iXi) + GPXi(iGPXi)*(Xi(iXi+1)-Xi(iXi)) )/2;
                        eta = ( Eta(iEta+1)+Eta(iEta) + GPEta(iGPEta)*(Eta(iEta+1)-Eta(iEta)) )/2;

                        %% 1iii.2ii. Find the knot span
                        xiSpan = findKnotSpan(xi,Xi,nxi);
                        etaSpan = findKnotSpan(eta,Eta,neta);

                        %% 1iii.2iii. Compute the NURBS basis functions for the patch
                        noDeriv = 1;
                        dR = computeIGABasisFunctionsAndDerivativesForSurface...
                            (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,noDeriv);

                        %% 1iii.2iv. Compute the physical location of the Gauss Point on patch
                        P = computeCartesianCoordinatesOfAPointOnBSplineSurface...
                            (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,dR(:,1));

                        %% 1iii.2v. Compute the parameter location of the Gauss Point on patch in patchReference
                        if isReferenceBSplinePatchProvided
                            if ~exist('xi0','var')
                                xi0 = xi0Init;
                            end
                            if ~exist('eta0','var')
                                eta0 = eta0Init;
                            end 
                            isConvergent = false;
                            [xiRef,etaRef,Projected,flagNR,~] = computeNearestPointProjectionOnBSplineSurface...
                                    (P,pRef,XiRef,qRef,EtaRef,CPRef,isNURBSRef,xi0,eta0,propNewtonRapshon);
                            if norm(Projected - P) < propError.tolClose || flagNR
                                isConvergent = true;
                            else
                                xi0 = Xi(1);
                                eta0 = Eta(1);
                            end
                            if ~isConvergent
                                while eta0 < Eta(end)
                                    while xi0 < Xi(end)
                                        [xiRef,etaRef,Projected,flagNR,~] = computeNearestPointProjectionOnBSplineSurface...
                                            (P,pRef,XiRef,qRef,EtaRef,CPRef,isNURBSRef,xi0,eta0,propNewtonRapshon);
                                        if norm(Projected - P) < propError.tolClose || flagNR
                                            isConvergent = true;
                                            break;
                                        end
                                        xi0 = xi0 + dxi0;
                                    end
                                    if isConvergent
                                        break;
                                    end
                                    xi0 = Xi(1);
                                    eta0 = eta0 + deta0;
                                end
                            end
                            if ~isConvergent
                                xi0 = xi0Init;
                                eta0 = eta0Init;
        %                         warning('Closest point projection on NURBS surface failed');
                                continue;
                            else
        %                         hold on;
        %                         plot3(P(1,1),P(2,1),P(3,1),'ws--','MarkerEdgeColor','g','MarkerFaceColor','g');
                            end
                        end

                        %% 1iii.2vi. Find the knot span indices for the reference patch
                        if isReferenceBSplinePatchProvided
                            xiSpanRef = findKnotSpan(xiRef,XiRef,nxiRef);
                            etaSpanRef = findKnotSpan(etaRef,EtaRef,netaRef);
                        end
                        
                        %% 1iii.2vii. Compute the NURBS basis functions for the reference patch
                        if isReferenceBSplinePatchProvided
                            noDeriv = 0;
                            RRef = computeIGABasisFunctionsAndDerivativesForSurface...
                                (xiSpanRef,pRef,xiRef,XiRef,etaSpanRef,qRef,etaRef,EtaRef,CPRef,isNURBSRef,noDeriv);
                        end

                        %% 1iii.2viii. Compute the base vectors at the Gauss points for the patch

                        % Compute the in-plane base vectors
                        noDrv = 0;
                        [A1,A2] = ...
                            computeBaseVectorsAndDerivativesForBSplineSurface...
                            (xiSpan,p,etaSpan,q,CP,noDrv,dR);

                        % Compute the not normalized surface normal vector
                        A3Tilde = cross(A1(:,1),A2(:,1));

                        % Compute the Jacobian of the transformation from the physical
                        % to the parameter space
                        detJPhys2Param = norm(A3Tilde);

                        %% 1iii.2ix. Compute the Cartesian coordinates of the Gauss points

                        % For patch :
                        % ___________

                        XPatch = computeCartesianCoordinatesOfAPointOnBSplineSurface...
                            (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,dR(:,1));

                        % For patchReference :
                        % ____________________

                        if isReferenceBSplinePatchProvided
                            XReference = computeCartesianCoordinatesOfAPointOnBSplineSurface...
                                (xiSpanRef,pRef,xiRef,XiRef,etaSpanRef,qRef,etaRef,EtaRef,CPRef,RRef);
                        elseif isReferenceSolutionProvided
                            XReference = computeReferenceSolution...
                                (XPatch,propReferenceSolution);
                        end

                        %% 1iii.2x. Compute the element area at the GP
                        elementAreaOnGP = detJPhys2Param*detJxizeta*GWXi(iGPXi)*GWEta(iGPEta);

                        %% 1iii.2xi. Compute the L2-norms

                        % Compute the L2-norms of the differences
                        diffGeoDomainL2(iPatches,1) = diffGeoDomainL2(iPatches,1) + ...
                            norm(XPatch - XReference)^2*elementAreaOnGP;

                        % Compute the L2-norms of the reference fields
                        geoRefDomainL2(iPatches,1) = geoRefDomainL2(iPatches,1) + ...
                            norm(XReference)^2*elementAreaOnGP;
                    end
                end
            end
        end
    end
end

%% 2. Loop over all the interfaces in the multipatch geometry
for iConnections = 1:noConnections
    %% 2i. Get the patch IDs
    
    % Patch I :
    % _________
    
    IDI = connections.xiEtaCoup(iConnections,1);
    
    % Patch J :
    % _________
    
    IDJ = connections.xiEtaCoup(iConnections,2);
    
    %% 2ii. Get the patch parameters
    
    % Patch I :
    % _________
    
    pI = BSplinePatches{IDI}.p;
    qI = BSplinePatches{IDI}.q;
    XiI = BSplinePatches{IDI}.Xi;
    EtaI = BSplinePatches{IDI}.Eta;
    CPI = BSplinePatches{IDI}.CP;
    isNURBSI = BSplinePatches{IDI}.isNURBS;
    nxiI = length(CPI(:,1,1));
    netaI = length(CPI(1,:,1));
    xicoupI = connections.xiEtaCoup(iConnections,3:4);
    etacoupI = connections.xiEtaCoup(iConnections,5:6);    
    
	% Patch 2 :
    % _________
    
    pJ = BSplinePatches{IDJ}.p;
    qJ = BSplinePatches{IDJ}.q;
    XiJ = BSplinePatches{IDJ}.Xi;
    EtaJ = BSplinePatches{IDJ}.Eta;
    CPJ = BSplinePatches{IDJ}.CP;
    isNURBS2 = BSplinePatches{IDJ}.isNURBS;
    nxiJ = length(CPJ(:,1,1));
    netaJ = length(CPJ(1,:,1));
    xicoupJ = connections.xiEtaCoup(iConnections,7:8);
    etacoupJ = connections.xiEtaCoup(iConnections,9:10);
    
    %% 2iii. Check if the patches have the same or opposite orientation
    haveSameOrientation = findSubdomainInterfaceOrientation...
        (pI,XiI,qI,EtaI,CPI,isNURBSI,xicoupI,etacoupI,pJ,XiJ,qJ,EtaJ,CPJ,isNURBS2,xicoupJ,etacoupJ);
    
    %% 2iv. Select an integration scheme
    if strcmp(propInt.type,'default')
        p = max(pI,pJ);
        q = max(qI,qJ);
        pol = max(p,q);
        noGP = ceil((pol + 1)/2);
    elseif strcmp(propInt.type,'user')
        noGP = propInt.noGP;
    else
        error('Define type in the propInt array for the interface integration');
    end
	[GP,GW] = getGaussPointsAndWeightsOverUnitDomain(noGP);
    
    %% 2v. Get the running and the fixed parameters on the patch interface and the coupling region

    % For patch I :
    % _____________

    if etacoupI(1) == etacoupI(2)
        % Coupled region in xi-direction
        couplingRegionI = xicoupI;

        % Find the correct spans for the coupled region
        spanStartI = findKnotSpan(couplingRegionI(1),XiI,nxiI);
        spanEndI = findKnotSpan(couplingRegionI(2),XiI,nxiI)+1;

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
        spanEndI = findKnotSpan(couplingRegionI(2),EtaI,netaI)+1;   

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
    
    %% 2vi. Compute the merged knot vector from both patches over the interface
    knotVector = unique(mergesorted(couplingRegionOnKnotVectorI,couplingRegionOnKnotVectorJ));
    
    %% 2vii. Loop over the elements on the coupling surface
    for i = 1:length(knotVector)-1
        if knotVector(i) ~= knotVector(i+1)
            %% 2vii.1. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
            detJxizeta = (knotVector(i+1)-knotVector(i))/2;
            
            %% 2vii.2. Loop over the Gauss points
            for iGP = 1:noGP
                %% 2vii.2i. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
                xiEta = ((1-GP(iGP))*knotVector(i)+(1+GP(iGP))*knotVector(i+1))/2;
                
                %% 2vii.2ii. Compute the NURBS basis functions
            
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
                    if ~haveSameOrientation
                        xiJ = XiJ(length(XiJ)) - xiJ;
                    end
                    xiSpanJ = findKnotSpan(xiJ,XiJ,nxiJ);
                else
                    etaJ = xiEta;
                    if ~haveSameOrientation
                        etaJ = EtaJ(length(EtaJ)) - etaJ;
                    end
                    etaSpanJ = findKnotSpan(etaJ,EtaJ,netaJ);
                end
                RJ = computeIGABasisFunctionsAndDerivativesForSurface...
                    (xiSpanJ,pJ,xiJ,XiJ,etaSpanJ,qJ,etaJ,EtaJ,CPJ,isNURBS2,0);
                
                %% 2vii.2iii. Compute the edge length on the Gauss Point
                [A1,A2] = computeBaseVectorsAndDerivativesForBSplineSurface...
                    (xiSpanI,pI,etaSpanI,qI,CPI,0,dRI);
                if isOnXiI
                    dL = length(A1);
                else
                    dL = length(A2);
                end
                
                %% 2vii.2iv. Compute the physical location of the Gauss point
                P = computeCartesianCoordinatesOfAPointOnBSplineSurface...
                    (xiSpanI,pI,xiI,XiI,etaSpanI,qI,etaI,EtaI,CPI,dRI(:,1));

                %% 2vii.2v. Compute the parameter location of the Gauss Point on patch in patchReference
                if isReferenceBSplinePatchProvided
                    if ~exist('xi0','var')
                        xi0 = xi0Init;
                    end
                    if ~exist('eta0','var')
                        eta0 = eta0Init;
                    end 
                    isConvergent = false;
                    [xiRef,etaRef,Projected,flagNR,~] = computeNearestPointProjectionOnBSplineSurface...
                            (P,pRef,XiRef,qRef,EtaRef,CPRef,isNURBSRef,xi0,eta0,propNewtonRapshon);
                    if norm(Projected - P) < propError.tolClose || flagNR
                        isConvergent = true;
                    else
                        xi0 = Xi(1);
                        eta0 = Eta(1);
                    end
                    if ~isConvergent
                        while eta0 < Eta(end)
                            while xi0 < Xi(end)
                                [xiRef,etaRef,Projected,flagNR,~] = computeNearestPointProjectionOnBSplineSurface...
                                    (P,pRef,XiRef,qRef,EtaRef,CPRef,isNURBSRef,xi0,eta0,propNewtonRapshon);
                                if norm(Projected - P) < propError.tolClose || flagNR
                                    isConvergent = true;
                                    break;
                                end
                                xi0 = xi0 + dxi0;
                            end
                            if isConvergent
                                break;
                            end
                            xi0 = Xi(1);
                            eta0 = eta0 + deta0;
                        end
                    end
                    if ~isConvergent
                        xi0 = xi0Init;
                        eta0 = eta0Init;
    %                         warning('Closest point projection on NURBS surface failed');
                        continue;
                    else
    %                         hold on;
    %                         plot3(P(1,1),P(2,1),P(3,1),'ws--','MarkerEdgeColor','g','MarkerFaceColor','g');
                    end
                end
                
                %% 2iii.2vi. Find the knot span indices for the reference patch
                if isReferenceBSplinePatchProvided
                    xiSpanRef = findKnotSpan(xiRef,XiRef,nxiRef);
                    etaSpanRef = findKnotSpan(etaRef,EtaRef,netaRef);
                end

                %% 2iii.2vii. Compute the NURBS basis functions for the reference patch
                if isReferenceBSplinePatchProvided
                    noDeriv = 0;
                    RRef = computeIGABasisFunctionsAndDerivativesForSurface...
                        (xiSpanRef,pRef,xiRef,XiRef,etaSpanRef,qRef,etaRef,EtaRef,CPRef,isNURBSRef,noDeriv);
                end
                
                %% 2vii.2viii. Compute the Cartesian coordinates of the Gauss points

                % For patch I :
                % _____________

                XPatchI = computeCartesianCoordinatesOfAPointOnBSplineSurface...
                    (xiSpanI,pI,xiI,XiI,etaSpanI,qI,etaI,EtaI,CPI,dRI(:,1));
                
                % For patch J :
                % _____________

                XPatchJ = computeCartesianCoordinatesOfAPointOnBSplineSurface...
                    (xiSpanJ,pJ,xiJ,XiJ,etaSpanJ,qJ,etaJ,EtaJ,CPJ,RJ(:,1));
                
                % For patchReference :
                % ____________________
                
                if isReferenceBSplinePatchProvided
                    XReference = computeCartesianCoordinatesOfAPointOnBSplineSurface...
                        (xiSpanRef,pRef,xiRef,XiRef,etaSpanRef,qRef,etaRef,EtaRef,CPRef,RRef);
                elseif isReferenceSolutionProvided
                    XReference = computeReferenceSolution...
                        (XPatch,propReferenceSolution);
                end
                
                %% 2vii.2ix. Compute the element length at the GP
                elementLengthOnGP = detJxizeta*dL*GW(iGP);
                
                %% 2vii.2xi. Compute the L2-norms
                
                % Compute the L2-norm of the differences
                diffGeoInterfaceL2(iConnections,1) = diffGeoInterfaceL2(iConnections,1) + ...
                    norm(XPatchI - XPatchJ)^2*elementLengthOnGP;
                
                % Compute the L2-norm of the exact geometry
                geoInterfaceL2(iConnections,1) = geoInterfaceL2(iConnections,1) + ...
                    norm(XReference)^2*elementLengthOnGP;
            end
        end
    end
end

%% 3. Get the relative values
if sqrt(geoRefDomainL2) > tol_null
    relErrorL2GeoDomain = sqrt(diffGeoDomainL2)./sqrt(geoRefDomainL2);
else
    relErrorL2GeoDomain = sqrt(diffGeoDomainL2);
end
if sqrt(geoInterfaceL2) > tol_null
    relErrorL2GeoInterface = sqrt(diffGeoInterfaceL2)./sqrt(geoInterfaceL2);
else
    relErrorL2GeoInterface = sqrt(diffGeoInterfaceL2);
end

%% 4. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;
    
    fprintf('Computation of the broken norm constituents took %.2d seconds \n\n',computationalTime);
    fprintf('______________________Linear Analysis Ended__________________________\n');
    fprintf('#####################################################################\n\n\n');
    fprintf('\n');
end

end
