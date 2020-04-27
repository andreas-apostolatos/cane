function [relErrorL2GeoDomain, relErrorL2GeoInterface, index] = ...
    computeDomainAndIterfaceErrorInL2AgainstFEMMembraneFormFinding ...
    (BSplinePatches, connections, mesh, propInt, graph, outMsg)
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
% for each interface between the patches with respect to the solution when
% classical Finite elements are used.
%
%                  Input :
%         BSplinePatches : Array of B-Spline surface patches
%            connections : Define the connection between the patches:
%                            .No : Number of connections
%                     .xiEtaCoup : [patchID1 patchID2 xi12 eta12 xi21 eta21
%                                   ...      ...    ...   ...  ...   ...]
%                   mesh : The form-found Finite element mesh consising of,
%                             .nodes : Array containing the nodal
%                                      coordinates and their global
%                                      numbering
%                          .elements : Array containing the elements
%                propInt : On the numerical integration over the domain :
%                            .type : 'default' or 'user'
%                         .noGPsXi : Number of Gauss points in xi-direction
%                        .noGPsEta : Number of Gauss points in eta-
%                                    direction
%                  graph : Structure containing the information on the
%                          figures,
%                            .index : The id of the figure to be produced
%                                     next
%                 outMsg : Enables outputting information onto the command
%                          window when chosen as 'outputEnabled'
%
%                 Output :
%    relErrorL2GeoDomain : Array containing the relative geometry error in 
%                          the domain for each patch
% relErrorL2GeoInterface : Array containing the geometrical jump across
%                          each interface
%                  index : The id of the next figure
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
%                1iii.2v. Find the closest to the Finite Element mesh node to the Gauss point
%
%               1iii.2vi. Find all elements containing the closest to the Gauss point node
%
%              1iii.2vii. Loop over all elements containing the closest to the Gauss point node and try to perform a projection
%
%             1iii.2viii. Compute the geometrical location of Gauss point in the Finite Element mesh
%
%               1iii.2ix. Compute the base vectors at the Gauss points for the patch
%
%                1iii.2x. Compute the element area at the GP
%
%               1iii.2xi. Compute the L2-norm of the errors
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
% 3. Loop over all patches, compute and visualize the error
% ->
%    3i. Get the patch parameters
%
%   3ii. Get the increment in the parameter space
%
%  3iii. Loop over all the parametric coordinates in eta-direction
%  ->
%        3iii.1. Find the span in the eta-direction
%
%        3iii.2. Initialize coordinate in xi-direction
%
%        3iii.3. Initialize counter in xi-direction
%
%        3iii.4. Loop over all the parametric coordinates in xi-direction
%        ->
%                3iii.4i. Find the span in xi-direction
%
%               3iii.4ii. Compute the IGA basis functions and possibly their derivatives
%
%              3iii.4iii. Compute and store the Cartesian coordinates of the parametric location
%
%               3iii.4iv. Find the closest to the Finite Element mesh node to the Gauss point
%
%                3iii.4v. Find all elements containing the closest to the Gauss point node
%
%               3iii.4vi. Loop over all elements containing the closest to the Gauss point node and try to perform a projection
%
%              3iii.4vii. Compute the geometrical location of Gauss point in the Finite Element mesh
%
%             3iii.4viii. Compute and store the geometrical error in the parametric location
%
%               3iii.4ix. Update the counter in xi-direction
%
%                3iii.4x. Update the parametric coordinate in xi-direction
%        <-
%
%        3iii.5. Update the counter in eta-direction
%
%        3iii.6. Update the parametric coordinate in eta-direction
%  <-
%
%   3iv. Visualize the geometrical error for the patch
% <-
%
% 4. Get the relative values
%
% 5. Appendix
%
%% function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('_____________________________________________________________________\n');
    fprintf('#####################################################################\n');
    fprintf('Computation of the constituents of the broken norm corresponding to \n');
    fprintf('the form-finding analysis of a multipatch membrane has been initiated \n');
    fprintf('_____________________________________________________________________\n');
    fprintf('\n');
    
    % start measuring computational time
    tic;
end

%% 0. Read input

% Check input
if ~isfield(mesh,'nodes')
    error('Variable mesh.nodes undefined');
end
if ~isfield(mesh,'elements')
    error('Variable mesh.elements undefined');
end

% Define a tolerance for considering the denominator equal to zero
tol_null = 1e-12;

% Define a tolerance for jumping out of the parameter space
tol = 1e-10;

% Number of B-Spline patches
noPatches = length(BSplinePatches);

% number of connections
noConnections = connections.No;

% Initialize auxiliary arrays for the figure production
isFigure = false;
if ~ischar(graph)
    if isfield(graph,'index')
        isFigure = true;
    end
    xiGrid = 49;
    etaGrid = 49;
end

% Initialize output values
diffGeoDomainL2 = zeros(noPatches,1);
geoRefDomainL2 = zeros(noPatches,1);
diffGeoInterfaceL2 = zeros(connections.No,1);

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
    for iEta = q + 1:meta - q -1
        for iXi = p + 1:mxi - p - 1
            if Xi(iXi+1) ~= Xi(iXi) && Eta(iEta + 1) ~= Eta(iEta)
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
                        XPatch = computeCartesianCoordinatesOfAPointOnBSplineSurface...
                            (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,dR(:,1));
                        
                        %% 1iii.2v. Find the closest to the Finite Element mesh node to the Gauss point
                        indexNode = dsearchn(mesh.nodes(:,2:4),XPatch');
                        idNode = mesh.nodes(indexNode,1);
                        
                        %% 1iii.2vi. Find all elements containing the closest to the Gauss point node
                        [indexElement,~] = find(mesh.elements == idNode);
                        
                        %% 1iii.2vii. Loop over all elements containing the closest to the Gauss point node and try to perform a projection
                        for iElm = 1:length(indexElement)
                            P1 = mesh.nodes(mesh.elements(indexElement(iElm,1),1),2:4)';
                            P2 = mesh.nodes(mesh.elements(indexElement(iElm,1),2),2:4)';
                            P3 = mesh.nodes(mesh.elements(indexElement(iElm,1),3),2:4)';
                            [~,N,isInside] = computePointProjectionAndBasisFunctionsOnLinearTriangle...
                                (XPatch,P1,P2,P3);
                            if isInside
                                break;
                            end
                        end
                        if ~isInside
                            continue;
                        end
                        
                        %% 1iii.2viii. Compute the geometrical location of Gauss point in the Finite Element mesh
                        XReference = [P1 P2 P3]*N;
                        
                        %% 1iii.2ix. Compute the base vectors at the Gauss points for the patch

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
                        
                        %% 1iii.2x. Compute the element area at the GP
                        elementAreaOnGP = detJPhys2Param*detJxizeta*GWXi(iGPXi)*GWEta(iGPEta);
                        
                        %% 1iii.2xi. Compute the L2-norm of the errors
                        
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
                
                %% 2vii.2iv. Compute the Cartesian coordinates of the Gauss points

                % For patch I :
                % _____________

                XPatchI = computeCartesianCoordinatesOfAPointOnBSplineSurface...
                    (xiSpanI,pI,xiI,XiI,etaSpanI,qI,etaI,EtaI,CPI,dRI(:,1));
                
                % For patch J :
                % _____________

                XPatchJ = computeCartesianCoordinatesOfAPointOnBSplineSurface...
                    (xiSpanJ,pJ,xiJ,XiJ,etaSpanJ,qJ,etaJ,EtaJ,CPJ,RJ(:,1));
                
                %% 2vii.2vi. Compute the element length at the GP
                elementLengthOnGP = detJxizeta*dL*GW(iGP);
                
                %% 2vii.2vii. Compute the L2-norm of the difference
                diffGeoInterfaceL2(iConnections,1) = diffGeoInterfaceL2(iConnections,1) + ...
                    norm(XPatchI - XPatchJ)^2*elementLengthOnGP;
            end
        end
    end
end

%% 3. Loop over all patches, compute and visualize the error
if isFigure
    figure(graph.index)
    for iPatches = 1:noPatches
        %% 3i. Get the patch parameters
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
        
        %% 3ii. Get the increment in the parameter space
        dxi = (Xi(end) - Xi(1))/(etaGrid - 1);
        deta = (Eta(end) - Eta(1))/(xiGrid - 1);
        
        %% 3iii. Initialize visualization arrays
        clear P errDomain;
        P = zeros(xiGrid,etaGrid,3);
        errDomain = zeros(xiGrid,etaGrid);
        counterEta = 1;
        eta = Eta(1);
        
        %% 3iii. Loop over all the parametric coordinates in eta-direction
        while eta <= Eta(meta) + tol
            %% 3iii.1. Find the span in the eta-direction
            etaSpan = findKnotSpan(eta,Eta,neta);

            %% 3iii.2. Initialize coordinate in xi-direction
            xi = Xi(1);

            %% 3iii.3. Initialize counter in xi-direction
            counterXi = 1;
            
            %% 3iii.4. Loop over all the parametric coordinates in xi-direction
            while xi <= Xi(mxi)+tol
                %% 3iii.4i. Find the span in xi-direction
                xiSpan = findKnotSpan(xi,Xi,nxi);

                %% 3iii.4ii. Compute the IGA basis functions and possibly their derivatives
                if strcmp(graph.resultant,'displacement')
                    noDrvBasis = 0;
                elseif strcmp(graph.resultant,'strain') || strcmp(graph.resultant,'force')
                    noDrvBasis = 1;
                end
                dR = computeIGABasisFunctionsAndDerivativesForSurface...
                    (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,noDrvBasis);
                
                %% 3iii.4iii. Compute and store the Cartesian coordinates of the parametric location
                 XPatch = computeCartesianCoordinatesOfAPointOnBSplineSurface...
                    (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,dR(:,1));
                P(counterXi,counterEta,:) = XPatch;
                
                %% 3iii.4iv. Find the closest to the Finite Element mesh node to the Gauss point
                indexNode = dsearchn(mesh.nodes(:,2:4),XPatch');
                idNode = mesh.nodes(indexNode,1);

                %% 3iii.4v. Find all elements containing the closest to the Gauss point node
                [indexElement,~] = find(mesh.elements == idNode);

                %% 3iii.4vi. Loop over all elements containing the closest to the Gauss point node and try to perform a projection
                for iElm = 1:length(indexElement)
                    P1 = mesh.nodes(mesh.elements(indexElement(iElm,1),1),2:4)';
                    P2 = mesh.nodes(mesh.elements(indexElement(iElm,1),2),2:4)';
                    P3 = mesh.nodes(mesh.elements(indexElement(iElm,1),3),2:4)';
                    [~,N,isInside] = computePointProjectionAndBasisFunctionsOnLinearTriangle...
                        (XPatch,P1,P2,P3);
                    if isInside
                        break;
                    end
                end
                
                %% 3iii.4vii. Compute the geometrical location of Gauss point in the Finite Element mesh
                if isInside
                    XReference = [P1 P2 P3]*N;
                else
                    XReference = mesh.nodes(indexNode,2:4)';
                end
                
                %% 3iii.4viii. Compute and store the geometrical error in the parametric location
                if isInside
                    errDomain(counterXi,counterEta) = norm(XPatch - XReference);
                else
                    errDomain(counterXi,counterEta) = 0;
                end
                
                %% 3iii.4ix. Update the counter in xi-direction
                counterXi = counterXi + 1;
                
                %% 3iii.4x. Update the parametric coordinate in xi-direction
                xi = xi + dxi;
            end
            
            %% 3iii.5. Update the counter in eta-direction
            counterEta = counterEta + 1;
            
            %% 3iii.6. Update the parametric coordinate in eta-direction
            eta = eta + deta;
        end
        
        %% 3iv. Visualize the geometrical error for the patch
        surf(P(:,:,1),P(:,:,2),P(:,:,3),errDomain);
        hold on;
        plot_knotsForBSplineSurfaceOnCartesianSpace(p,q,Xi,Eta,CP,isNURBS,false,xiGrid,etaGrid);
    end
    hold off;
    shading interp;
%     colormap('jet');
%     colorbar;
    % invert default colormap => red = negativ, blue = positive
    % COL = colormap;
    % invCOL(:,1) = COL(:,3);
    % invCOL(:,2) = COL(:,2);
    % invCOL(:,3) = COL(:,1);
    % colormap(invCOL);
    % make colormap symmetric
    % colim = caxis;
    % caxis([-max(abs(colim)) max(abs(colim))]);
    axis equal;
    index = graph.index + 1;
else
    index = 'undefined'; 
end
                        
%% 4. Get the relative values
if sqrt(geoRefDomainL2) > tol_null
    relErrorL2GeoDomain = sqrt(diffGeoDomainL2)./sqrt(geoRefDomainL2);
else
    relErrorL2GeoDomain = sqrt(diffGeoDomainL2);
end
relErrorL2GeoInterface = sqrt(diffGeoInterfaceL2);

%% 5. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;
    
    fprintf('Computation of the broken norm constituents took %.2d seconds \n\n',computationalTime);
    fprintf('______________________Linear Analysis Ended__________________________\n');
    fprintf('#####################################################################\n\n\n');
    fprintf('\n');
end

end
