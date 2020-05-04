function KWeakDBCPenaltyMembrane = ...
    computeWeakDBCMtxPenaltyIGAMembrane ...
    (BSplinePatch, connections, numDOFs, propCoupling)
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
% boundary conditions weakly using the Penalty method.
%
%               Input :
%        BSplinePatch : The B-Spline patch array containing:
%                      .p,q : The polynomial orders of the B-Spline surface
%                             in both parametric directions
%                   .Xi,Eta : The knot vectors in both parametric
%                             directions
%                       .CP : The set of control points and weights
%                  .isNURBS : Flag on whether the basis is a NURBS or a
%                             B-Spline
%               .parameters : Technical parameters for the structure
%                  .homDOFs : The global numbering of the DOFs where 
%                             homogeneous Dirichlet boundary conditions are
%                             applied
%                .inhomDOFs : The global numbering of the DOFs where
%                             inhomogeneous Dirichlet boundary conditions
%                             are applied
%          .valuesInhomDOFs : The values on the DOFs corresponding to the
%                             application of inhomogeneous Dirichlet
%                             boundary conditions
%                      .NBC : Structure containing information on the 
%                             application of the Neumann boundary 
%                             conditions
%                                       .noCnd : Number of Neumann boundary 
%                                                conditions
%                             .xiLoadExtension : Cell array {.noCnd} 
%                                                containing the load 
%                                                extensions in the xi-
%                                                direction
%                            .etaLoadExtension : Cell array {.noCnd} 
%                                                containing the load 
%                                                extensions in the eta-
%                                                direction
%                               .loadAmplitude : Array (1,.noCnd) 
%                                                containing the load 
%                                                amplitudes
%                               .loadDirection : Array (1,.noCnd) 
%                                                containing the load 
%                                                directions
%                              .computeLoadVct : Cell array {.noCnd} 
%                                                containing the function 
%                                                name for the computation 
%                                                of the load vector
%                              .isConservative : Array (1,.noCnd) of flags 
%                                                indicating whether the 
%                                                load is conservative or 
%                                                not
%                  .weakDBC : Structure containing information on the
%                             application of weak boundary conditions:
%                                       .alpha : The penalty parameter
%                                 .xiExtension : xi-extension of the
%                                                boundary where the weak
%                                                boundary conditions are
%                                                applied
%                                .etaExtension : eta-extension of the
%                                                boundary where the weak
%                                                boundary conditions are
%                                                applied
%         connections : Dummy variable for this function
%             numDOFs : Dummy variable for this function
%        propCoupling : Dummy variable for this function
%
% Function layout :
%
% 0. Read input
% 
% 1. Loop over all the conditions
% ->
%    1i. Get the extensions of the Dirichlet boundary where weak boundary conditions using the penalty method are to be applied
%
%   1ii. Get the running and the fixed parameters on the patch boundary where the weak boundary conditions are applied
%
%  1iii. Get the parameter space of the application of the weak boundary conditions using the knot vector information
%
%   1iv. Issue Gauss Point coordinates and weights
%
%    1v. Loop over the elements on the parameter space where the weak boundary conditions are applied
%    ->
%        1v.1. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
%
%        1v.2. Loop over the Gauss points
%        ->
%              1v.2i. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
%
%             1v.2ii. Compute the NURBS basis functions
%
%            1v.2iii. Create the element freedom table
%
%             1v.2iv. Compute the covariant base vectors
%
%              1v.2v. Compute the B-operator matrix for the displacement field
%
%             1v.2vi. Compute the determinant of the Jacobian of the transformation from the physical to the parent domain on the GP
%
%            1v.2vii. Compute the element length at the GP
%
%           1v.2viii. Compute the element stiffness contribution matrix for the weak application of the Dirichlet boundary conditions and add it to the global matrix
%        <-
%    <-
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
numCPs_xi = length(CP(:, 1, 1));
numCPs_eta = length(CP(1, :, 1));

% Number of local DOFs
noNodesLoc = (p + 1)*(q + 1);
noDOFsLoc = 3*noNodesLoc;

% Number of DOFs
numDOFs = 3*numCPs_xi*numCPs_eta;

% Initialize the element freedom table
EFT = zeros(1,noDOFsLoc);

% Initialize auxiliary arrays
RMtx = zeros(3,noDOFsLoc);

% Initialize the output arrays
KWeakDBCPenaltyMembrane = zeros(numDOFs,numDOFs);

%% 1. Loop over all the conditions
for counterCnd = 1:BSplinePatch.weakDBC.noCnd
    %% 1i. Get the extensions of the Dirichlet boundary where weak boundary conditions using the penalty method are to be applied
    xiExtension = BSplinePatch.weakDBC.xiExtension{counterCnd};
    etaExtension = BSplinePatch.weakDBC.etaExtension{counterCnd};
    isWeakDBCOverPoint = false;
    
    %% 1ii. Get the running and the fixed parameters on the patch boundary where the weak boundary conditions are applied
    if etaExtension(1) == etaExtension(2) && xiExtension(1) ~= xiExtension(2)
        % Coupled region in xi-direction
        couplingRegion = xiExtension;

        % Find the correct spans for the coupled region
        spanStart = findKnotSpan(couplingRegion(1),Xi,numCPs_xi);
        spanEnd = findKnotSpan(couplingRegion(2),Xi,numCPs_xi) + 1;

        % Find the corresponding to the coupled region knot span
        weakDBCRgionOnKnotVector = Xi(spanStart:spanEnd);

        % Fixed parameter on the parametric net
        eta = etaExtension(1);

        % Find the span where xiEta it lies in
        etaSpan = findKnotSpan(eta,Eta,numCPs_eta);

        % Flag on whether the coupling line is over xi
        isOnXi = true;
    elseif xiExtension(1) == xiExtension(2) && etaExtension(1) ~= etaExtension(2)
        % Coupled region in eta-direction
        couplingRegion = etaExtension;

        % Find the correct spans for the coupled region
        spanStart = findKnotSpan(couplingRegion(1),Eta,numCPs_eta);   
        spanEnd = findKnotSpan(couplingRegion(2),Eta,numCPs_eta) + 1;

        % Find the corresponding to the coupled region knot span
        weakDBCRgionOnKnotVector = Eta(spanStart:spanEnd);

        % Fixed parameter on the parametric net
        xi = xiExtension(1);

        % Find the span where uv it lies in
        xiSpan = findKnotSpan(xi,Xi,numCPs_xi);

        % Flag on whether the coupling line is over eta
        isOnXi = false;
    elseif etaExtension(1) == etaExtension(2) && xiExtension(1) == xiExtension(2)
        % Initialize flag on whether the application of weak Dirichlet
        % boundary conditions is over a point or over a line
        isWeakDBCOverPoint = true;
        
        % Fixed parameters on the parametric net
        xi = xiExtension(1);
        eta = etaExtension(1);
        
        % Find the correct spans for the point
        xiSpan = findKnotSpan(xi,Xi,numCPs_xi);
        etaSpan = findKnotSpan(eta,Eta,numCPs_eta);
        
        % Find the corresponding to the coupled region knot span
        weakDBCRgionOnKnotVector = zeros(2,1);
    else
        error('Boundary over which weak Dirichlet boundary conditions are to be imposed wrongly defined');
    end

    %% 1iii. Get the parameter space of the application of the weak boundary conditions using the knot vector information
    if ~isWeakDBCOverPoint
        weakDBCRgionOnKnotVector = unique(weakDBCRgionOnKnotVector);
    end
    
    %% 1iv. Issue Gauss Point coordinates and weights
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
%     GP = fliplr(GP);
%     GW = fliplr(GW);

    %% 1v. Loop over the elements on the parameter space where the weak boundary conditions are applied
    for i = 1:length(weakDBCRgionOnKnotVector)-1
        if (weakDBCRgionOnKnotVector(i) ~= weakDBCRgionOnKnotVector(i+1)) || isWeakDBCOverPoint
            %% 1v.1. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
            if ~isWeakDBCOverPoint
                detJxizeta = (weakDBCRgionOnKnotVector(i+1)-weakDBCRgionOnKnotVector(i))/2;
            end

            %% 1v.2. Loop over the Gauss points
            for j = 1:noGPs
                %% 1v.2i. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
                if ~isWeakDBCOverPoint
                    xiEta = ((1-GP(j))*weakDBCRgionOnKnotVector(i) + (1+GP(j))*weakDBCRgionOnKnotVector(i+1))/2;
                end

                %% 1v.2ii. Compute the NURBS basis functions
                if ~isWeakDBCOverPoint
                    if isOnXi
                        xi = xiEta;
                        xiSpan = findKnotSpan(xi,Xi,numCPs_xi);
                    else
                        eta = xiEta;
                        etaSpan = findKnotSpan(eta,Eta,numCPs_eta);
                    end
                end
                dR = computeIGABasisFunctionsAndDerivativesForSurface...
                    (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,1);

                %% 1v.2iii. Create the element freedom table

                % Initialize of the counter
                r = 1;

                % Relation global-local DoFs
                for cpj = etaSpan-q:etaSpan
                    for cpi = xiSpan-p:xiSpan
                        EFT(r)   = DOFNumbering(cpi,cpj,1);
                        EFT(r+1) = DOFNumbering(cpi,cpj,2);
                        EFT(r+2) = DOFNumbering(cpi,cpj,3);

                        % update counter
                        r = r + 3;
                    end
                end

                %% 1v.2iv. Compute the covariant base vectors
                [A1,A2] = computeBaseVectorsAndDerivativesForBSplineSurface...
                    (xiSpan,p,etaSpan,q,CP,0,dR);

                %% 1v.2v. Compute the B-operator matrix for the displacement field

                % initialize counter
                k = 0;

                % Loop over all the non-zero contributions at the span
                % under study
                for c = 0:q
                    for b = 0:p
                        % Update counter
                        k = k + 1;

                        % Matrix containing the basis functions
                        RMtx(1,3*k-2) = dR(k,1);
                        RMtx(2,3*k-1) = dR(k,1);
                        RMtx(3,3*k) = dR(k,1);
                    end
                end

                %% 1v.2vi. Compute the determinant of the Jacobian of the transformation from the physical to the parent domain on the GP
                if ~isWeakDBCOverPoint
                    if isOnXi
                        detJxxi = norm(A1(:,1));
                    else
                        detJxxi = norm(A2(:,1));
                    end
                end

                %% 1v.2vii. Compute the element length at the GP
                if ~isWeakDBCOverPoint
                    elementLengthOnGP = detJxxi*detJxizeta*GW(j);
                else
                    elementLengthOnGP = 1;
                end

                %% 1v.2viii. Compute the element stiffness contribution matrix for the weak application of the Dirichlet boundary conditions and add it to the global matrix
                KWeakDBCPenaltyMembrane(EFT,EFT) = KWeakDBCPenaltyMembrane(EFT,EFT) + ...
                    BSplinePatch.weakDBC.alpha*(RMtx'*RMtx)*elementLengthOnGP;
            end
        end
    end
end

end
