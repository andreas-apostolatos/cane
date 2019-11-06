function [errorL2_weakDBC,errorL2_diffInterfaceDisp,errorL2_diffInterfaceTrac] = ...
    computeErrorsIGAMembraneMultipatchesNLinearTransient...
    (BSplinePatches,connections,propTransientAnalysis,dHatHistory,int)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the jump of interface quantities as well as domain quantities 
% when a reference field is provided. The function returns also the errors
% along the interface where weakly imposed Dirichlet conditions are assumed
%
%                    Input :
%           BSplinePatches : Structure containing all the information 
%                            regarding the connections between the 
%                            multipatches
%              connections : Define the connection between the patches:
%                            .No : Number of connections
%                     .xiEtaCoup : [patchID1 patchID2 xi12 eta12 xi21 eta21
%                                 ...      ...    ...   ...  ...   ...]
%    propTransientAnalysis : On the transient analysis
%                           .TStart : Start time of the simulation
%                             .TEnd : End time of the simulation
%                      .noTimeSteps : Number of time steps of the
%                                     simulation
%                               .dt : Time step of the simulation
%              dHatHistory : The history of the displacement field for the
%                            multipatch geometry
%                      int : Quadrature rule for the computation of the
%                            interface jumps and relative domain errors :
%                               .type : 'default' or 'user'
%                               .noGP : Number of Gauss points for the
%                                       error integration if 'user'-type of
%                                       quadrature is selected
%
%                   Output :
%          errorL2_weakDBC :  Array [noTimeSteps,noCnd,3] which contains 
%                             the displacement error at the Dirichlet 
%                             boundaries where the conditions are weakly 
%                             imposed in the L2-norm as well as the patch 
%                             the condition belongs to and the numbering of 
%                             the condition within that patch
% errorL2_diffInterfaceDisp : Array [noTimeSteps,noConnections,3] which
%                             contains the displacement error along the
%                             interfaces of the multipatch geometry in the
%                             L2-norm as well as the numbering of the
%                             patches involved
% errorL2_diffInterfaceTrac : Array [noTimeSteps,noConnections,3] which
%                             contains the traction error along the
%                             interfaces of the multipatch geometry in the
%                             L2-norm as well as the numbering of the
%                             patches involved
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the time instance of the transient simulation
% ->
%    1i. Get actual time
%
%   1ii. Loop over the interface connections in the multipatch geometry
%   ->
%        1ii.1. Get the patch IDs
%
%        1ii.2. Get the coupling boundaries for both patches
%
%        1ii.3. Determine the interface orientation
%
%        1ii.4. Recover the B-Spline patch parameters
%
%        1ii.5. Get the displacement vector for each patch and compute the displaced Control Points
%
%        1ii.6. Get the running and the fixed parameters on the patch interface and the coupling region
%
%        1ii.7. Compute the merged knot vector from both patches over the interface
%
%        1ii.8. Choose an interface integration rule
%
%        1ii.9. Loop over the elements on the coupling interface
%        ->
%               1ii.9i. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
%
%               1ii.9ii. Loop over the Gauss points
%               ->
%                        1ii.9ii.1. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
%
%                        1ii.9ii.2. Compute the NURBS basis functions
%
%                        1ii.9ii.3. Create the element freedom tables
%
%                        1ii.9ii.4. Compute the covariant base vectors of the reference configuration
%
%                        1ii.9ii.5. Compute the normal to the boundary vector and transform it into the contavariant basis
%
%                        1ii.9ii.6. Compute the covariant base vectors of the current configuration
%
%                        1ii.9ii.7. Compute the covariant metric coefficients of the current configuration
%
%                        1ii.9ii.8. Compute the contravariant bases
%
%                        1ii.9ii.9. Compute the local Cartesian bases
%
%                        1ii.9ii.10. Compute the transformation matrix from the contravariant basis to the local Cartesian one
%
%                        1ii.9ii.11. Compute the transformation matrix from the local Cartesian bases to the covariant one
%
%                        1ii.9ii.12. Compute the Green-Lagrange strains in the contravariant bases
%
%                        1ii.9ii.13. Transform the Green-Lagrange strains in the local Cartesian bases
%
%                        1ii.9ii.14. Compute the prestress values on the local Cartesian coordinate systems
%
%                        1ii.9ii.15. Compute the 2nd Piola-kirchhoff stresses in the local Cartesian systems
%
%                        1ii.9ii.16. Transform the 2nd Piola-kirchhoff stresses in the covariant systems
%
%                        1ii.9ii.17. Compute the stress components
%
%                        1ii.9ii.18. Compute the traction vectors
%
%                        1ii.9ii.19. Compute the basis functions matrix and their derivatives and the determinant of the Jacobian to the transformation from the physical space (x-y) to the NURBS parameter space (xi-eta)
%
%                        1ii.9ii.20. Get the displacement field at the Gauss point for each patch
%
%                        1ii.9ii.21. Compute the determinant of the Jacobian of the transformation from the physical to the parent domain on the GP
%
%                        1ii.9ii.22. Compute the element length at the GP
%
%                        1ii.9ii.23. Compute the difference of the displacements and tractions from each patch at the Gauss point and add the contribution
%               <-
%        <-
%    <-
%
%  1iii. Initialize the counter of the boundaries where weakly Dirichlet boundary conditions are applied
%
%   1iv. Loop over the patches of the multipatch geometry
%   ->
%        1iv.1. Retrieve the patch information
%      
%        1iv.2. Get the displacement vector for each patch at the current time step
%
%        1iv.3. Loop over the boundaries of the patch where weakly Dirichlet boundary conditions are applied
%        ->
%                1iv.3i. Assign the numbering of the corresponding patch and the corresponding Dirichlet boundary condition
%
%               1iv.3ii. Get the boundary of the patch where weakly Dirichlet boundary conditions are applied
%
%              1iv.3iii. Check along which parametric direction the Dirichlet boundary conditions are weakly applied
%
%               1iv.3iv. Choose an boundary integration rule
%
%                1iv.3v. Loop over the elements on the coupling interface
%                ->
%                        1iv.3v.1. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
%
%                        1iv.3v.2. Loop over the Gauss points
%                        ->
%                                   1iv.3v.2i. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
%
%                                  1iv.3v.2ii. Compute the NURBS basis functions
%
%                                 1iv.3v.2iii. Create the element freedom tables
%
%                                  1iv.3v.2iv. Compute the covariant base vectors of the reference configuration
%
%                                   1iv.3v.2v. Compute the basis functions matrix and their derivatives and the determinant of the Jacobian to the transformation from the physical space (x-y) to the NURBS parameter space (xi-eta)
%
%                                  1iv.3v.2vi. Get the displacement field at the Gauss point
%
%                                 1iv.3v.2vii. Get the prescribed displacement field at the Gauss point
%
%                                1iv.3v.2viii. Compute the determinant of the Jacobian of the transformation from the physical to the parent domain on the GP
%
%                                  1iv.3v.2ix. Compute the element length at the GP
%
%                                   1iv.3v.2x. Compute the difference of the displacements and tractions from each patch at the Gauss point and add the contribution
%                        <-
%                <-
%
%               1iv.3vi. Update the counter of the conditions
%        <-
%   <-
% <-
%
% 3. Get the square root of the interface quantities
%
%% Function main body

%% 0. Read input

% Number of time steps
noTimeSteps = length(dHatHistory(1,:));

% Define tolerances
epsArea = 1e-15;

% Get the number of patches
noPatches = length(BSplinePatches);

% Check input
for iPatches = 1:noPatches
    if ~isfield(BSplinePatches{iPatches},'EFTPatches')
        error('Variable BSplinePatches{%d}.EFTPatches is undefined. Ensure that BSplinePatches is written out by a solver',iPatches);
    end
end
if noTimeSteps ~= propTransientAnalysis.noTimeSteps + 2
    error('The number of time steps and the transient results do not match');
end

% Check the number of the Dirichlet boundaries where weakly imposed
% conditions are assumed
noCnd = 0;
for iPatches = 1:noPatches
    for iCnd = 1:BSplinePatches{iPatches}.weakDBC.noCnd
        noCnd = noCnd + 1;
    end
end

% Initialize output arrays
errorL2_weakDBC = zeros(noTimeSteps,noCnd,3);
errorL2_diffInterfaceDisp = zeros(noTimeSteps,connections.No,3);
errorL2_diffInterfaceTrac = zeros(noTimeSteps,connections.No,3);

%% 1. Loop over all the time instance of the transient simulation
for iTimeSteps = 1:noTimeSteps
    %% 1i. Get actual time
    t = propTransientAnalysis.TStart + (iTimeSteps - 1)*propTransientAnalysis.dt;
    
    %% 1ii. Loop over the interface connections in the multipatch geometry
    for iConnections = 1:connections.No
        %% 1ii.1. Get the patch IDs

        % Patch I :
        % _________

        IDI = connections.xiEtaCoup(iConnections,1);
        errorL2_diffInterfaceDisp(iTimeSteps,iConnections,2) = IDI;
        errorL2_diffInterfaceTrac(iTimeSteps,iConnections,2) = IDI;

        % Patch J :
        % _________

        IDJ = connections.xiEtaCoup(iConnections,2);
        errorL2_diffInterfaceDisp(iTimeSteps,iConnections,3) = IDJ;
        errorL2_diffInterfaceTrac(iTimeSteps,iConnections,3) = IDJ;

        %% 1ii.2. Get the coupling boundaries for both patches

        % Patch I :
        % _________

        BSplinePatches{IDI}.xicoup = connections.xiEtaCoup(iConnections,3:4);
        BSplinePatches{IDI}.etacoup = connections.xiEtaCoup(iConnections,5:6);

        % Patch J :
        % _________

        BSplinePatches{IDJ}.xicoup = connections.xiEtaCoup(iConnections,7:8);
        BSplinePatches{IDJ}.etacoup = connections.xiEtaCoup(iConnections,9:10);

        %% 1ii.3. Determine the interface orientation
        haveSameOrientation = findSubdomainInterfaceOrientation...
            (BSplinePatches{IDI}.p,BSplinePatches{IDI}.Xi,BSplinePatches{IDI}.q,BSplinePatches{IDI}.Eta,BSplinePatches{IDI}.CP,BSplinePatches{IDI}.isNURBS,BSplinePatches{IDI}.xicoup,BSplinePatches{IDI}.etacoup,...
            BSplinePatches{IDJ}.p,BSplinePatches{IDJ}.Xi,BSplinePatches{IDJ}.q,BSplinePatches{IDJ}.Eta,BSplinePatches{IDJ}.CP,BSplinePatches{IDJ}.isNURBS,BSplinePatches{IDJ}.xicoup,BSplinePatches{IDJ}.etacoup);

        %% 1ii.4. Recover the B-Spline patch parameters

        % Patch I :
        % _________

        pI = BSplinePatches{IDI}.p;
        qI = BSplinePatches{IDI}.q;
        XiI = BSplinePatches{IDI}.Xi;
        EtaI = BSplinePatches{IDI}.Eta;
        CPI = BSplinePatches{IDI}.CP;
        isNURBSI = BSplinePatches{IDI}.isNURBS;
        parametersI = BSplinePatches{IDI}.parameters;
        thicknessI = parametersI.t;
        materialMtxVoigtI = parametersI.E*thicknessI/(1-parametersI.nue^2)*...
            [1               parametersI.nue 0
             parametersI.nue 1               0
             0               0               (1-parametersI.nue)/2];
        prestressI = parametersI.prestress;
        xicoupI = BSplinePatches{IDI}.xicoup;
        etacoupI = BSplinePatches{IDI}.etacoup;
        nxiI = length(CPI(:,1,1));
        netaI = length(CPI(1,:,1));
        DOFNumberingI = BSplinePatches{IDI}.DOFNumbering;
        noNodesLocI = (pI + 1)*(qI + 1);
        noDOFsLocI = 3*noNodesLocI;
        BDisplacementsGCI = zeros(3,noDOFsLocI);
        EFTI = zeros(1,noDOFsLocI);

        % Patch J :
        % _________

        pJ = BSplinePatches{IDJ}.p;
        qJ = BSplinePatches{IDJ}.q;
        XiJ = BSplinePatches{IDJ}.Xi;
        EtaJ = BSplinePatches{IDJ}.Eta;
        CPJ = BSplinePatches{IDJ}.CP;
        isNURBSJ = BSplinePatches{IDJ}.isNURBS;
        parametersJ = BSplinePatches{IDJ}.parameters;
        thicknessJ = parametersJ.t;
        materialMtxVoigtJ = parametersJ.E*thicknessJ/(1-parametersJ.nue^2)*...
            [1               parametersJ.nue 0
             parametersJ.nue 1               0
             0               0               (1-parametersJ.nue)/2];
        prestressJ = parametersJ.prestress;
        xicoupJ = BSplinePatches{IDJ}.xicoup;
        etacoupJ = BSplinePatches{IDJ}.etacoup;
        nxiJ = length(CPJ(:,1,1));
        netaJ = length(CPJ(1,:,1));
        DOFNumberingJ = BSplinePatches{IDJ}.DOFNumbering;
        noNodesLocJ = (pJ + 1)*(qJ + 1);
        noDOFsLocJ = 3*noNodesLocJ;
        BDisplacementsGCJ = zeros(3,noDOFsLocJ);
        EFTJ = zeros(1,noDOFsLocJ);

        %% 1ii.5. Get the displacement vector for each patch and compute the displaced Control Points

        % For patch I :
        % _____________

        dHatI = dHatHistory(BSplinePatches{IDI}.EFTPatches,iTimeSteps);
        CPdI = computeDisplacedControlPointsForIGAKirchhoffLoveShell...
                        (CPI,dHatI);
                    
        % For patch J :
        % _____________
        
        dHatJ = dHatHistory(BSplinePatches{IDJ}.EFTPatches,iTimeSteps);
        CPdJ = computeDisplacedControlPointsForIGAKirchhoffLoveShell...
                        (CPJ,dHatJ);
        
        %% 1ii.6. Get the running and the fixed parameters on the patch interface and the coupling region

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

        %% 1ii.7. Compute the merged knot vector from both patches over the interface

        % Merge the two knot vectors into one for integration purposes:
        couplingRegionOnKnotVector = mergesorted(couplingRegionOnKnotVectorI,couplingRegionOnKnotVectorJ);

        % Delete double entries
        couplingRegionOnKnotVector = unique(couplingRegionOnKnotVector);

        %% 1ii.8. Choose an interface integration rule
        if strcmp(int.type,'default')
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
        elseif strcmp(int.type,'user')
            noGPs = int.nGP;
        end
        [GP,GW] = getGaussPointsAndWeightsOverUnitDomain(noGPs);
        
        %% 1ii.9. Loop over the elements on the coupling interface
        for iElmts = 1:length(couplingRegionOnKnotVector)-1
            if couplingRegionOnKnotVector(iElmts) ~= couplingRegionOnKnotVector(iElmts+1)
                %% 1ii.9i. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
                detJxizeta = (couplingRegionOnKnotVector(iElmts+1) - couplingRegionOnKnotVector(iElmts))/2;

                %% 1ii.9ii. Loop over the Gauss points
                for iGPs = 1:noGPs
                    %% 1ii.9ii.1. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
                    xiEta = ((1-GP(iGPs))*couplingRegionOnKnotVector(iElmts) + ...
                        (1+GP(iGPs))*couplingRegionOnKnotVector(iElmts+1))/2;
                    
                    %% 1ii.9ii.2. Compute the NURBS basis functions
                    
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
                    dRJ = computeIGABasisFunctionsAndDerivativesForSurface...
                        (xiSpanJ,pJ,xiJ,XiJ,etaSpanJ,qJ,etaJ,EtaJ,CPJ,isNURBSJ,1);
                    
                    %% 1ii.9ii.3. Create the element freedom tables
            
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
                    
                    %% 1ii.9ii.4. Compute the covariant base vectors of the reference configuration
                    
                    % For patch I :
                    % _____________

                    [A1I,A2I] = computeBaseVectorsAndDerivativesForBSplineSurface...
                        (xiSpanI,pI,etaSpanI,qI,CPI,0,dRI);
                    A3I = cross(A1I,A2I)/norm(cross(A1I,A2I));

                    % For patch J :
                    % _____________

                    [A1J,A2J] = computeBaseVectorsAndDerivativesForBSplineSurface...
                        (xiSpanJ,pJ,etaSpanJ,qJ,CPJ,0,dRJ);
                    A3J = cross(A1J,A2J)/norm(cross(A1J,A2J));
                    
                    %% 1ii.9ii.5. Compute the normal to the boundary vector and transform it into the contavariant basis
                    
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
                    
                    %% 1ii.9ii.6. Compute the covariant base vectors of the current configuration
            
                    % For patch I :
                    % _____________

                    [a1I,a2I] = computeBaseVectorsAndDerivativesForBSplineSurface...
                        (xiSpanI,pI,etaSpanI,qI,CPdI,0,dRI);

                    % For patch J :
                    % _____________

                    [a1J,a2J] = computeBaseVectorsAndDerivativesForBSplineSurface...
                        (xiSpanJ,pJ,etaSpanJ,qJ,CPdJ,0,dRJ);
                    
                    %% 1ii.9ii.7. Compute the covariant metric coefficients of the current configuration
            
                    % For patch I :
                    % _____________

                    aabCovariantI = [a1I a2I]'*[a1I a2I];

                    % For patch J :
                    % _____________

                    aabCovariantJ = [a1J a2J]'*[a1J a2J];
                    
                    %% 1ii.9ii.8. Compute the contravariant bases
            
                    % For patch I :
                    % _____________

                    AabCovariantI = [A1I A2I]'*[A1I A2I];
                    if det(AabCovariantI) < epsArea
                        continue;
                    end
                    AContravariantI = AabCovariantI\[A1I A2I]';
                    AContravariantI = AContravariantI';

                    % For patch J :
                    % _____________

                    AabCovariantJ = [A1J A2J]'*[A1J A2J];
                    if det(AabCovariantJ) < epsArea
                        continue;
                    end
                    AContravariantJ = AabCovariantJ\[A1J A2J]';
                    AContravariantJ = AContravariantJ';
                    
                    %% 1ii.9ii.9. Compute the local Cartesian bases
            
                    % For patch I :
                    % _____________

                    eLCI = computeLocalCartesianBasis4BSplineSurface...
                        ([A1I A2I],AContravariantI);

                    % For patch J :
                    % _____________

                    eLCJ = computeLocalCartesianBasis4BSplineSurface...
                        ([A1J A2J],AContravariantJ);
                    
                    %% 1ii.9ii.10. Compute the transformation matrix from the contravariant basis to the local Cartesian one
            
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
                    
                    %% 1ii.9ii.11. Compute the transformation matrix from the local Cartesian bases to the covariant one
            
                    % For patch I :
                    % _____________

                    TFromLCToCovI = computeTFromLocalCartesian2CovariantBasis4BSplineSurface...
                        (eLCI,AContravariantI);

                    % For patch I :
                    % _____________

                    TFromLCToCovJ = computeTFromLocalCartesian2CovariantBasis4BSplineSurface...
                        (eLCJ,AContravariantJ);
                    
                    %% 1ii.9ii.12. Compute the Green-Lagrange strains in the contravariant bases
            
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
                             
                    %% 1ii.9ii.13. Transform the Green-Lagrange strains in the local Cartesian bases
            
                    % For patch I :
                    % _____________

                    EpsilonLCI = TFromContraToLC4VoigtStrainI*EpsilonContraI;

                    % For patch I :
                    % _____________

                    EpsilonLCJ = TFromContraToLC4VoigtStrainJ*EpsilonContraJ;
                    
                    %% 1ii.9ii.14. Compute the prestress values on the local Cartesian coordinate systems
            
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
                    
                    %% 1ii.9ii.15. Compute the 2nd Piola-kirchhoff stresses in the local Cartesian systems
            
                    % For patch I :
                    % _____________

                    NLCI = thicknessI*pTildeI +  materialMtxVoigtI*EpsilonLCI;

                    % For patch J :
                    % _____________

                    NLCJ = thicknessJ*pTildeJ +  materialMtxVoigtJ*EpsilonLCJ;
                    
                    %% 1ii.9ii.16. Transform the 2nd Piola-kirchhoff stresses in the covariant systems
            
                    % For patch I :
                    % _____________

                    NCovariantI = TFromLCToCovI*NLCI;

                    % For patch J :
                    % _____________

                    NCovariantJ = TFromLCToCovJ*NLCJ;
                    
                    %% 1ii.9ii.17. Compute the stress components
            
                    % For patch I :
                    % _____________

                    PalphaBetaI = [NCovariantI(1,1) NCovariantI(3,1)
                                   NCovariantI(3,1) NCovariantI(2,1)];

                    % For patch J :
                    % _____________

                    PalphaBetaJ = [NCovariantJ(1,1) NCovariantJ(3,1)
                                   NCovariantJ(3,1) NCovariantJ(2,1)];

                    %% 1ii.9ii.18. Compute the traction vectors
            
                    % For patch I :
                    % _____________

                    tractionVctI = [a1I a2I]*PalphaBetaI*uContravariantI;

                    % For patch J :
                    % _____________

                    tractionVctJ = [a1J a2J]*PalphaBetaJ*uContravariantJ;
                    
                    %% 1ii.9ii.19. Compute the basis functions matrix and their derivatives and the determinant of the Jacobian to the transformation from the physical space (x-y) to the NURBS parameter space (xi-eta)
                    
                    % For patch I :
                    % _____________

                    % initialize counter
                    kI = 0;

                    % Loop over all the non-zero contributions at the span
                    % under study
                    for c = 0:qI
                        for b = 0:pI
                            % Update counter
                            kI = kI + 1;

                            % Matrix containing the basis functions
                            BDisplacementsGCI(1,3*kI-2) = dRI(kI,1);
                            BDisplacementsGCI(2,3*kI-1) = dRI(kI,1);
                            BDisplacementsGCI(3,3*kI) = dRI(kI,1);
                        end
                    end

                    % For patch J :
                    % _____________

                    % initialize counter
                    kJ = 0;

                    % Loop over all the non-zero contributions at the span
                    % under study
                    for c = 0:qJ
                        for b = 0:pJ
                            % Update counter
                            kJ = kJ + 1;

                            % Matrix containing the basis functions
                            BDisplacementsGCJ(1,3*kJ-2) = dRJ(kJ,1);
                            BDisplacementsGCJ(2,3*kJ-1) = dRJ(kJ,1);
                            BDisplacementsGCJ(3,3*kJ) = dRJ(kJ,1);
                        end
                    end
                    
                    %% 1ii.9ii.20. Get the displacement field at the Gauss point for each patch
                    
                    % For patch I :
                    % _____________
                    
                    dispVctI = BDisplacementsGCI*dHatI(EFTI);
                    
                    % For patch J :
                    % _____________
                    
                    dispVctJ = BDisplacementsGCJ*dHatJ(EFTJ);
                    
                    %% 1ii.9ii.21. Compute the determinant of the Jacobian of the transformation from the physical to the parent domain on the GP
                    if isOnXiI
                        detJxxi = norm(A1I(:,1));
                    else
                        detJxxi = norm(A2I(:,1));
                    end
                    
                    %% 1ii.9ii.22. Compute the element length at the GP
                    elementLengthOnGP = detJxxi*detJxizeta*GW(iGPs);
            
                    %% 1ii.9ii.23. Compute the difference of the displacements and tractions from each patch at the Gauss point and add the contribution
                    errorL2_diffInterfaceTrac(iTimeSteps,iConnections,1) = ...
                        errorL2_diffInterfaceTrac(iTimeSteps,iConnections,1) + ...
                        norm(tractionVctI + tractionVctJ)*elementLengthOnGP;
                    errorL2_diffInterfaceDisp(iTimeSteps,iConnections,1) = ...
                        errorL2_diffInterfaceDisp(iTimeSteps,iConnections,1) + ...
                        norm(dispVctI - dispVctJ)*elementLengthOnGP;
                end
            end
        end
    end
    
    %% 1iii. Initialize the counter of the boundaries where weakly Dirichlet boundary conditions are applied
    counterCnd = 1;
    
    %% 1iiv. Loop over the patches of the multipatch geometry
    for iPatches = 1:noPatches
        %% 1iv.1. Retrieve the patch information
        p = BSplinePatches{iPatches}.p;
        q = BSplinePatches{iPatches}.q;
        Xi = BSplinePatches{iPatches}.Xi;
        Eta = BSplinePatches{iPatches}.Eta;
        CP = BSplinePatches{iPatches}.CP;
        isNURBS = BSplinePatches{iPatches}.isNURBS;
        weakDBC = BSplinePatches{iPatches}.weakDBC;
        nxi = length(CP(:,1,1));
        neta = length(CP(1,:,1));
        DOFNumbering = BSplinePatches{iPatches}.DOFNumbering;
        noNodesLoc = (p + 1)*(q + 1);
        noDOFsLoc = 3*noNodesLoc;
        BDisplacementsGC = zeros(3,noDOFsLoc);
        EFT = zeros(1,noDOFsLoc);
        
        %% 1iv.2. Get the displacement vector for each patch at the current time step
        dHat = dHatHistory(BSplinePatches{iPatches}.EFTPatches,iTimeSteps);
        
        %% 1iv.3. Loop over the boundaries of the patch where weakly Dirichlet boundary conditions are applied
        for iCnd = 1:weakDBC.noCnd
            %% 1iv.3i. Assign the numbering of the corresponding patch and the corresponding Dirichlet boundary condition
            errorL2_weakDBC(iTimeSteps,counterCnd,2) = iPatches;
            errorL2_weakDBC(iTimeSteps,counterCnd,3) = iCnd;
            
            %% 1iv.3ii. Get the boundary of the patch where weakly Dirichlet boundary conditions are applied
            xiWeakDBC = weakDBC.xiExtension{iCnd};
            etaWeakDBC = weakDBC.etaExtension{iCnd};
            
            %% 1iv.3iii. Check along which parametric direction the Dirichlet boundary conditions are weakly applied
            if etaWeakDBC(1) == etaWeakDBC(2)
                % Coupled region in xi-direction
                weakDBCRegion = xiWeakDBC;

                % Find the correct spans for the coupled region
                spanStart = findKnotSpan(weakDBCRegion(1),Xi,nxi);
                spanEnd = findKnotSpan(weakDBCRegion(2),Xi,nxi) + 1;

                % Corresponding to the coupled region knot span
                weakDBCRegionOnKnotVector = Xi(spanStart:spanEnd);

                % Fixed parameter on the parametric net
                eta = etaWeakDBC(1);

                % Find the span where xiEta it lies in
                etaSpan = findKnotSpan(eta,Eta,neta);

                % Flag on whether the coupling line is over xi
                isOnXi = true;
            else
                % Coupled region in eta-direction
                weakDBCRegion = etaWeakDBC;

                % Find the correct spans for the coupled region
                spanStart = findKnotSpan(weakDBCRegion(1),Eta,neta);   
                spanEnd = findKnotSpan(weakDBCRegion(2),Eta,neta) + 1;   

                % Corresponding to the coupled region knot span
                weakDBCRegionOnKnotVector = Eta(spanStart:spanEnd);

                % Fixed parameter on the parametric net
                xi = xiWeakDBC(1);

                % Find the span where uv it lies in
                xiSpan = findKnotSpan(xi,Xi,nxi);

                % Flag on whether the coupling line is over eta
                isOnXi = false;
            end
            
            %% 1iv.3iv. Choose an boundary integration rule
            if strcmp(int.type,'default')
                if isOnXi
                    noGPs = p + 1;
                else
                    noGPs = q + 1;
                end
            elseif strcmp(int.type,'user')
                noGPs = int.noGPs;
            end
            [GP,GW] = getGaussPointsAndWeightsOverUnitDomain(noGPs);
            
            %% 1iv.3v. Loop over the elements on the coupling interface
            for iElmts = 1:length(weakDBCRegionOnKnotVector)-1
                if weakDBCRegionOnKnotVector(iElmts) ~= weakDBCRegionOnKnotVector(iElmts + 1)
                    %% 1iv.3v.1. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
                    detJxizeta = (weakDBCRegionOnKnotVector(iElmts + 1) - weakDBCRegionOnKnotVector(iElmts))/2;

                    %% 1iv.3v.2. Loop over the Gauss points
                    for iGPs = 1:noGPs
                        %% 1iv.3v.2i. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
                        xiEta = ((1-GP(iGPs))*weakDBCRegionOnKnotVector(iElmts) + ...
                            (1+GP(iGPs))*weakDBCRegionOnKnotVector(iElmts + 1))/2;

                        %% 1iv.3v.2ii. Compute the NURBS basis functions
                        if isOnXi
                            xi = xiEta;
                            xiSpan = findKnotSpan(xi,Xi,nxi);
                        else
                            eta = xiEta;
                            etaSpan = findKnotSpan(eta,Eta,neta);
                        end
                        dR = computeIGABasisFunctionsAndDerivativesForSurface...
                            (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,1);
                        
                        %% 1iv.3v.2iii. Compute the Cartesian coordinates of the Gauss Point
                        if isfield(weakDBC,'imposedMotion')
                            if isa(weakDBC.imposedMotion{iCnd},'function_handle')
                                X = computeCartesianCoordinatesOfAPointOnBSplineSurface...
                                    (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,dR(:,1));
                            end
                        end

                        %% 1iv.3v.2iii. Create the element freedom tables
                        
                        % Initialize of the counter
                        r = 1;

                        % Relation global-local DoFs
                        for cpj = etaSpan - q:etaSpan
                            for cpi = xiSpan - p:xiSpan
                                EFT(r) = DOFNumbering(cpi,cpj,1);
                                EFT(r + 1) = DOFNumbering(cpi,cpj,2);
                                EFT(r + 2) = DOFNumbering(cpi,cpj,3);

                                % update counter
                                r = r + 3;
                            end
                        end

                        %% 1iv.3v.2iv. Compute the covariant base vectors of the reference configuration
                        [A1,A2] = computeBaseVectorsAndDerivativesForBSplineSurface...
                            (xiSpan,p,etaSpan,q,CP,0,dR);
                        
                        %% 1iv.3v.2v. Compute the basis functions matrix and their derivatives and the determinant of the Jacobian to the transformation from the physical space (x-y) to the NURBS parameter space (xi-eta)

                        % initialize counter
                        k = 0;

                        % Loop over all the non-zero contributions at the span
                        % under study
                        for c = 0:q
                            for b = 0:p
                                % Update counter
                                k = k + 1;

                                % Matrix containing the basis functions
                                BDisplacementsGC(1,3*k - 2) = dR(k,1);
                                BDisplacementsGC(2,3*k - 1) = dR(k,1);
                                BDisplacementsGC(3,3*k) = dR(k,1);
                            end
                        end

                        %% 1iv.3v.2vi. Get the displacement field at the Gauss point
                        dispVct = BDisplacementsGC*dHat(EFT);
                        
                        %% 1iv.3v.2vii. Get the prescribed displacement field at the Gauss point
                        if isfield(weakDBC,'imposedMotion')
                            if isa(weakDBC.imposedMotion{iCnd},'function_handle')
                                g2 = weakDBC.imposedMotion{iCnd}(X(1,1),X(2,1),X(3,1),t);
                            else
                                g2 = zeros(3,1);
                            end
                        else
                            g2 = zeros(3,1);
                        end
                        
                        %% 1iv.3v.2viii. Compute the determinant of the Jacobian of the transformation from the physical to the parent domain on the GP
                        if isOnXi
                            detJxxi = norm(A1(:,1));
                        else
                            detJxxi = norm(A2(:,1));
                        end
                        
                        %% 1iv.3v.2ix. Compute the element length at the GP
                        elementLengthOnGP = detJxxi*detJxizeta*GW(iGPs);
                        
                        %% 1iv.3v.2x. Compute the difference of the displacements and tractions from each patch at the Gauss point and add the contribution
                        errorL2_weakDBC(iTimeSteps,counterCnd,1) = ...
                            errorL2_weakDBC(iTimeSteps,counterCnd,1) + ...
                            norm(dispVct - g2)*elementLengthOnGP;
                    end
                end
            end
            %% 1iv.3vi. Update the counter of the conditions
            counterCnd = counterCnd + 1;
        end
    end
end

%% 3. Get the square root of the interface quantities
for iTimeSteps = 1:noTimeSteps
    for iConnections = 1:connections.No
        errorL2_diffInterfaceDisp(iTimeSteps,iConnections,1) = ...
            sqrt(errorL2_diffInterfaceDisp(iTimeSteps,iConnections,1));
        errorL2_diffInterfaceTrac(iTimeSteps,iConnections,1) = ...
            sqrt(errorL2_diffInterfaceTrac(iTimeSteps,iConnections,1));
    end
    for iCnd = 1:noCnd
        errorL2_weakDBC(iTimeSteps,iCnd,1) = sqrt(errorL2_weakDBC(iTimeSteps,iCnd,1));
    end
end

end
