function [trCurveNumGP, trCurveGPs, trCurveGPWeights, ...
    trCurveGPTangents, trCurveGPJacobianProducts] = ...
    computeGaussPointDataAlongDirichletBoundary...
    (BSplinePatch,xiExtension,etaExtension,propIntDirichlet,tab,outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the Gauss Point data along a given Dirichlet boundary.
%
%                     Input :
%              BSplinePatch : Structure array containing information on the
%                             B-Spline,
%                                       .p,q : The polynomial orders of the 
%                                              B-Spline surface in both 
%                                              parametric directions
%                                    .Xi,Eta : The knot vectors in both 
%                                              parametric directions
%                                        .CP : The set of control points 
%                                              and weights
%                                   .isNURBS : Flag on whether the basis is 
%                                              a NURBS or a B-Spline
%               xiExtension : Extension along the xi-direction of the
%                             Dirichlet boundary
%              etaExtension : Extension along the eta-direction of the
%                             Dirichlet boundary
%          propIntDirichlet : Integration properties,
%                               .type : 'default', 'user'
%                              .noGPs : Number of Gauss Points
%                       tab : Tabulation when printing message in the 
%                             command window
%                    outMsg : Enables message outputting in the command 
%                             window if chosen as 'outputEnabled'
%
%                    Output :
%              trCurveNumGP : Total number of Gauss points along the 
%                             Dirichlet boundary
%                trCurveGPs : Array containing the parametric locations of 
%                             all GPs along the Dirichlet boundary
%          trCurveGPWeights : Array containing all the Gauss weights
%         trCurveGPTangents : Array containing all the coordinates of the 
%                             tangent vectors in the Cartesian space along 
%                             the Dirichlet boundary
% trCurveGPJacobianProducts : Array containing all the Jacobian products
%                             times the Gauss point weights at the Gauss
%                             point along the Dirichlet boundary
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
    fprintf([tab 'Creating Gauss Point data along Dirichlet boundary\n']);
end

%% 0. Read input

% Check the NURBS geometry input
if iscell(BSplinePatch)
    if length(BSplinePatch) > 1
        error('Multipatch NURBS surface is given as input to the computation of the stiffness matrix for a single patch NURBS surface');
    else
        BSplinePatch = BSplinePatch{1};
    end
end

% Reassign the analysis arrays
p = BSplinePatch.p;
q = BSplinePatch.q;
Xi = BSplinePatch.Xi;
Eta = BSplinePatch.Eta;
CP = BSplinePatch.CP;
isNURBS = BSplinePatch.isNURBS;
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));

% Tolerance for the determinant of the covariant metric coefficient tensor
scaleTolDet = 1e-2;

% Flag on whether the integration domain is a point
isWeakDBCOverPoint = false;

% Initialize counters
counterGPs = 1;

%% 1. Get the running and the fixed parameters on the patch boundary where the weak boundary conditions are applied
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

%% 2. Issue the Gauss quadrature
if ~isWeakDBCOverPoint
    if strcmp(propIntDirichlet.type,'default')
        if isOnXi
            pDegree = 2*p;
        else
            pDegree = 2*q;
        end
        noGPs = ceil((pDegree + 1)/2);
    elseif strcmp(propIntDirichlet.type,'user')
        noGPs = propIntDirichlet.noGPs;
    end
    [GP,GW] = getGaussPointsAndWeightsOverUnitDomain(noGPs);
else
    noGPs = 1;
end
if strcmp(outMsg,'outputEnabled')
    fprintf([tab '\t' 'Creating %d Gauss Points per knot span\n'],noGPs);
end

%% 3. Compute the total number of Gauss points and initialize output arrays
noGPsTotal = 0;
for i = 1:length(weakDBCRgionOnKnotVector) - 1
    if (weakDBCRgionOnKnotVector(i) ~= weakDBCRgionOnKnotVector(i+1)) || isWeakDBCOverPoint
        noGPsTotal = noGPsTotal + noGPs;
    end
end
if strcmp(outMsg,'outputEnabled')
    fprintf([tab '\t' 'Creating data for %d Gauss Points in total\n'],noGPsTotal);
end
trCurveNumGP = noGPsTotal;
trCurveGPs = zeros(1,2*noGPsTotal);
trCurveGPWeights = zeros(1,noGPsTotal);
trCurveGPTangents = zeros(1,3*noGPsTotal);
trCurveGPJacobianProducts = zeros(1,noGPsTotal);

%% 4. Loop over the elements on the parameter space where the weak boundary conditions are applied
for i = 1:length(weakDBCRgionOnKnotVector) - 1
     if (weakDBCRgionOnKnotVector(i) ~= weakDBCRgionOnKnotVector(i+1)) || isWeakDBCOverPoint
        %% 4i. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
        if ~isWeakDBCOverPoint
            detJxizeta = (weakDBCRgionOnKnotVector(i+1) - weakDBCRgionOnKnotVector(i))/2;
        end
        
        %% 4ii. Loop over the Gauss points
        for j = 1:noGPs
            %% 4ii.1. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
            if ~isWeakDBCOverPoint
                xiEta = ((1-GP(j))*weakDBCRgionOnKnotVector(i) + (1+GP(j))*weakDBCRgionOnKnotVector(i+1))/2;
            end
            
            %% 4ii.2. Compute the NURBS basis functions
            if ~isWeakDBCOverPoint
                if isOnXi
                    xi = xiEta;
                    xiSpan = findKnotSpan(xi,Xi,nxi);
                else
                    eta = xiEta;
                    etaSpan = findKnotSpan(eta,Eta,neta);
                end
            end
            trCurveGPs(1,2*counterGPs - 1) = xi;
            trCurveGPs(1,2*counterGPs) = eta;
            dR = computeIGABasisFunctionsAndDerivativesForSurface...
                    (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,1);
                
            %% 4ii.3. Compute the covariant base vectors of the reference configuration
            [A1,A2] = computeBaseVectorsAndDerivativesForBSplineSurface...
                (xiSpan,p,etaSpan,q,CP,0,dR);
            if isOnXi
                trCurveGPTangents(1,3*counterGPs - 2) = A1(1,1);
                trCurveGPTangents(1,3*counterGPs - 1) = A1(2,1);
                trCurveGPTangents(1,3*counterGPs) = A1(3,1);
            else
                trCurveGPTangents(1,3*counterGPs - 2) = A2(1,1);
                trCurveGPTangents(1,3*counterGPs - 1) = A2(2,1);
                trCurveGPTangents(1,3*counterGPs) = A2(3,1);
            end
            
            %% 4ii.4. Compute the determinant of the Jacobian of the transformation from the physical to the parent domain on the GP
            if ~isWeakDBCOverPoint
                if isOnXi
                    detJxxi = norm(A1(:,1));
                else
                    detJxxi = norm(A2(:,1));
                end
                if detJxxi < tolDet
                    warning('Base vector with almost zero length was deteceted');
                end
            end
            
            %% 4ii.5. Compute the element length at the GP
            if ~isWeakDBCOverPoint
                elementLengthOnGP = detJxxi*detJxizeta*GW(j);
            else
                elementLengthOnGP = 1;
            end
            trCurveGPWeights(1,counterGPs) = GW(j);
            trCurveGPJacobianProducts(1,counterGPs) = elementLengthOnGP;
            
            %% 4ii.6. Update the Gauss Point counter
            counterGPs = counterGPs + 1;
        end
    end
end

end
