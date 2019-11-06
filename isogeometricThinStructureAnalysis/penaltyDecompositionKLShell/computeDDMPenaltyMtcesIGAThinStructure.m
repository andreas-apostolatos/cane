function [KpDispI,KpRotI,CpDispI,CpRotI,KpDispJ,KpRotJ] = ...
    computeDDMPenaltyMtcesIGAThinStructure...
    (patchI,patchJ,alphaD,alphaR,haveSameOrientation,int)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
%  Returns the stiffness and the coupling matrices for the penalty
%  decomposition method applied to the multipatch Kirchhoff-Love shell. The
%  coupling matrices are only return for patch I cause the respective
%  contributions in patch J are symmetric.
%
%                  Input : 
%          patchI,patchJ : The B-Spline patches which are sharing an
%                          interface
%                 alphaD : Penalty factor for the displacement coupling
%                 alphaR : Penalty factor for the rotation coupling
%    haveSameOrientation : Flag on whether the couling surfaces of the 
%                          shells are oriented in the same direction over 
%                          the coupling interface
%                   int : On the interface quadrature :
%                           .type : 'default' or 'user'
%                          .noGPs : Number of Gauss Points
%
%                Output :
%        KpDispI,KpRotI : Displacement and rotation coupling contribution
%                         to the stiffness matrix for patch 1
%        KpDispJ,KpRotJ : Displacement and rotation coupling contribution
%                         to the stiffness matrix for patch 2
%        CpDispI,CpRotI : Displacement and rotation coupling contribution
%                         to the coupling matrix for patch 1
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
% 4. Loop over the elements on the coupling interface
% ->
%    4i. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
%
%  4iii. Loop over the Gauss points
%  ->
%        4iii.1. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
%
%        4iii.2. Compute the NURBS basis functions
%
%        4iii.3. Create the element freedom tables for both patches
%
%        4iii.4. Compute the covariant base vectors
%
%        4iii.5. Compute the surface normal vectors
%
%        4iii.6. Compute the derivatives of the surface normal vectors
%
%        4iii.7. Compute the normal to the boundary vector
%
%        4iii.8. Compute the normal to the boundary vector
%
%        4iii.9. Compute the covariant metric coefficients
%
%        4iii.10. Compute the contravariant base vectors
%
%        4iii.11. Compute the basis functions matrix and their derivatives and the determinant of the Jacobian to the transformation from the physical space (x-y) to the NURBS parameter space (xi-eta)
%
%        4iii.12. Transform the normal and the tangent vectors to the covariant bases
%
%        4iii.13. Compute the curvature coefficients
%
%        4iii.14. Compute the B-operator matrices for the rotations
%
%        4iii.15. Compute the determinant of the Jacobian of the transformation from the physical to the parent domain on the GP
%
%        4iii.16. Compute the element length at the GP
%
%        4iii.17. Compute the element stiffness matrix contributions on the GP and add them to the global matrices
%
%        4iii.18. Compute the element coupling matrices at the GP
%  <-
% <-
% 
%% Function main body

%% 0. Read input

% Initialize a tolerance for the geometrical map
tolDet = 1e-8;

% For patch I :
% _____________

% Reassign the analysis arrays
pI = patchI.p;
qI = patchI.q;
XiI = patchI.Xi;
EtaI = patchI.Eta;
CPI = patchI.CP;
isNURBSI = patchI.isNURBS;
xicoupI = patchI.xicoup;
etacoupI = patchI.etacoup;
nDOFsI = patchI.noDOFs;

% Get the DOF numbering
DOFNumberingI = patchI.DOFNumbering;

% Number of Control Points in xi-,eta- directions
nxiI = length(CPI(:,1,1));
netaI = length(CPI(1,:,1));

% Number of local DOFs
nNodesLocI = (pI+1)*(qI+1);
nDOFsLocI = 3*nNodesLocI;

% For patch J :
% _____________

pJ = patchJ.p;
qJ = patchJ.q;
XiJ = patchJ.Xi;
EtaJ = patchJ.Eta;
CPJ = patchJ.CP;
isNURBSJ = patchJ.isNURBS;
xicoupJ = patchJ.xicoup;
etacoupJ = patchJ.etacoup;
nDOFsJ = patchJ.noDOFs;

% Get the DOF numbering
DOFNumberingJ = patchJ.DOFNumbering;

% Number of Control Points in xi-,eta- directions
nxiJ = length(CPJ(:,1,1));
netaJ = length(CPJ(1,:,1));

% Number of local DOFs
nNodesLocJ = (pJ+1)*(qJ+1);
nDOFsLocJ = 3*nNodesLocJ;

% Initialize the element freedom tables
EFTI = zeros(1,nDOFsLocI);
EFTJ = zeros(1,nDOFsLocJ);

% Initialize auxiliary arrays
BDisplacementsGCI = zeros(3,nDOFsLocI);
BDisplacementsGCJ = zeros(3,nDOFsLocJ);
dRdxiI = zeros(3,nDOFsLocI);
dRdxiJ = zeros(3,nDOFsLocJ);
dRdetaI = zeros(3,nDOFsLocI);
dRdetaJ = zeros(3,nDOFsLocJ);

% Initialize the output arrays
KpDispI = zeros(nDOFsI,nDOFsI);
KpRotI = zeros(nDOFsI,nDOFsI);
KpDispJ = zeros(nDOFsJ,nDOFsJ);
KpRotJ = zeros(nDOFsJ,nDOFsJ);
CpDispI = zeros(nDOFsI,nDOFsJ);
CpRotI = zeros(nDOFsI,nDOFsJ);

%% 1. Get the running and the fixed parameters on the patch interface and the coupling region

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

% Merge the two knot vectors into one for integration purposes:
couplingRegionOnKnotVector = mergesorted(couplingRegionOnKnotVectorI,couplingRegionOnKnotVectorJ);

% Delete double entries
couplingRegionOnKnotVector = unique(couplingRegionOnKnotVector);

%% 3. Issue Gauss Point coordinates and weights
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
    noGPs = int.noGPs;
end
[GP,GW] = getGaussPointsAndWeightsOverUnitDomain(noGPs);
GP = fliplr(GP);
GW = fliplr(GW);

%% 4. Loop over the elements on the coupling interface
for i = 1:length(couplingRegionOnKnotVector)-1
    if couplingRegionOnKnotVector(i) ~= couplingRegionOnKnotVector(i+1)
        %% 4i. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
        detJxizeta = (couplingRegionOnKnotVector(i+1)-couplingRegionOnKnotVector(i))/2;

        %% 4iii. Loop over the Gauss points
        for j = 1:noGPs
            %% 4iii.1. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
            xiEta = ((1-GP(j))*couplingRegionOnKnotVector(i)+(1+GP(j))*couplingRegionOnKnotVector(i+1))/2;

            %% 4iii.2. Compute the NURBS basis functions
            
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
                (xiSpanI,pI,xiI,XiI,etaSpanI,qI,etaI,EtaI,CPI,isNURBSI,2);
            
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
                (xiSpanJ,pJ,xiJ,XiJ,etaSpanJ,qJ,etaJ,EtaJ,CPJ,isNURBSJ,2);
            
            %% 4iii.3. Create the element freedom tables
            
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
        
            %% 4iii.4. Compute the covariant base vectors
            
            % For patch I :
            % _____________
            
            [dA1I,dA2I] = computeBaseVectorsAndDerivativesForBSplineSurface...
                (xiSpanI,pI,etaSpanI,qI,CPI,1,dRI);
            
            % For patch J :
            % _____________
            
            [dA1J,dA2J] = computeBaseVectorsAndDerivativesForBSplineSurface...
                (xiSpanJ,pJ,etaSpanJ,qJ,CPJ,1,dRJ);
            
            %% 4iii.5. Compute the surface normal vectors
            
            % For patch I :
            % _____________
            
            A3TildeI = cross(dA1I(:,1),dA2I(:,1));
            A3I = A3TildeI/norm(A3TildeI);
            
            % For patch J :
            % _____________
            
            A3TildeJ = cross(dA1J(:,1),dA2J(:,1));
            A3J = A3TildeJ/norm(A3TildeJ);
            
            % Assign a flag on the surface normal orientation
            A3ITimesA3J = A3I'*A3J;
            if A3ITimesA3J < 0
                haveSurfaceNormalsSameOrientation = false;
            else
                haveSurfaceNormalsSameOrientation = true;
            end
            if norm(A3TildeI) < tolDet || norm(A3TildeJ) < tolDet
                continue;
            end
            
            %% 4iii.6. Compute the derivatives of the surface normal vectors
            
            % For patch I :
            % _____________
            
            [dA3I,~] = computeParametricDrvsSurfaceNormalOnBSplineSurface...
                ([dA1I(:,1) dA2I(:,1)],[dA1I(:,2) dA2I(:,2) dA1I(:,3)],...
                A3I,norm(A3TildeI));
            
            % For patch J :
            % _____________
            
            [dA3J,~] = computeParametricDrvsSurfaceNormalOnBSplineSurface...
                ([dA1J(:,1) dA2J(:,1)],[dA1J(:,2) dA2J(:,2) dA1J(:,3)],...
                A3J,norm(A3TildeJ));
            
            %% 4iii.7. Compute the normal to the boundary vector
            
            % For patch I :
            % _____________
            
            [nI,tI] = computeNormalAndTangentVectorsToBSplineBoundary...
                (xiI,XiI,etaI,EtaI,dA1I(:,1),dA2I(:,1),A3I,isOnXiI);
            
            % For patch J :
            % _____________
            
            [nJ,tJ] = computeNormalAndTangentVectorsToBSplineBoundary...
                (xiJ,XiJ,etaJ,EtaJ,dA1J(:,1),dA2J(:,1),A3J,isOnXiJ);
            
            %% 4iii.8. Compute the covariant metric coefficients
            
            % For patch I :
            % _____________
            
            AabCovI = [dA1I(:,1) dA2I(:,1)]'*[dA1I(:,1) dA2I(:,1)];
            
            % For patch J :
            % _____________
            
            AabCovJ = [dA1J(:,1) dA2J(:,1)]'*[dA1J(:,1) dA2J(:,1)];
            
            %% 4iii.9. Compute the contravariant base vectors
            
            % For patch I :
            % _____________
            
            AContravariantI = (AabCovI\[dA1I(:,1) dA2I(:,1)]')';
            
            % For patch J :
            % _____________
            
            AContravariantJ = (AabCovJ\[dA1J(:,1) dA2J(:,1)]')';
            
            %% 4iii.10. Compute the basis functions matrix and their derivatives and the determinant of the Jacobian to the transformation from the physical space (x-y) to the NURBS parameter space (xi-eta)
            
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

                    % Matrix containing the derivatives of the basis functions
                    % With respect to xi:
                    dRdxiI(1,3*kI-2) = dRI(kI,2);
                    dRdxiI(2,3*kI-1) = dRI(kI,2);
                    dRdxiI(3,3*kI) = dRI(kI,2);

                    % With respect to eta:
                    dRdetaI(1,3*kI-2) = dRI(kI,4);
                    dRdetaI(2,3*kI-1) = dRI(kI,4);
                    dRdetaI(3,3*kI) = dRI(kI,4);
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

                    % Matrix containing the derivatives of the basis functions
                    % With respect to xi:
                    dRdxiJ(1,3*kJ-2) = dRJ(kJ,2);
                    dRdxiJ(2,3*kJ-1) = dRJ(kJ,2);
                    dRdxiJ(3,3*kJ) = dRJ(kJ,2);

                    % With respect to eta:
                    dRdetaJ(1,3*kJ-2) = dRJ(kJ,4);
                    dRdetaJ(2,3*kJ-1) = dRJ(kJ,4);
                    dRdetaJ(3,3*kJ) = dRJ(kJ,4);
                end
            end
            
            %% 4iii.11. Transform the normal and the tangent vectors to the covariant bases
            
            % For patch I :
            % _____________
            
            nCovariantI = AContravariantI'*nI;
            tCovariantI = AContravariantI'*tI;
            
            % For patch J :
            % _____________
            
            nCovariantJ = AContravariantJ'*nJ;
            tCovariantJ = AContravariantJ'*tJ;
            
            %% 4iii.12. Compute the curvature coefficients
            
            % For patch I :
            % _____________
            
            BVI = [dA1I(:,2) dA2I(:,2) dA1I(:,3)]'*A3I;
            
            % For patch J :
            % _____________
            
            BVJ = [dA1J(:,2) dA2J(:,2) dA1J(:,3)]'*A3J;
            
            %% 4iii.13. Compute the B-operator matrices for the rotations
            
            % For patch I :
            % _____________
            
            [BtI,BnI,~,~] = computeBOperatorMatrix4RotationsIGAKirchhoffLoveShell...
                (BDisplacementsGCI,dRdxiI,dRdetaI,A3I,dA3I,...
                AContravariantI,BVI,nCovariantI,tCovariantI);
            BRotationsI = [BtI
                           BnI];
            
            % For patch J :
            % _____________
            
            [BtJ,BnJ,~,~] = computeBOperatorMatrix4RotationsIGAKirchhoffLoveShell...
                (BDisplacementsGCJ,dRdxiJ,dRdetaJ,A3J,dA3J,...
                AContravariantJ,BVJ,nCovariantJ,tCovariantJ);
            BRotationsJ = [BtJ
                           BnJ];
                       
            %% 4iii.14. Compute the determinant of the Jacobian of the transformation from the physical to the parent domain on the GP
            if isOnXiI
                detJxxi = norm(dA1I(:,1));
            else
                detJxxi = norm(dA2I(:,1));
            end
            
            %% 4iii.16. Compute the element length at the GP
            elementLengthOnGP = detJxxi*detJxizeta*GW(j);
            
            %% 4iii.17. Compute the element stiffness matrix contributions on the GP and add them to the global matrices
            
            % For patch I :
            % _____________
            
            % Compute the displacement stiffness matrix
            if ~strcmp(alphaD,'undefined')
                KpDispI(EFTI,EFTI) = KpDispI(EFTI,EFTI) + ...
                    alphaD*(BDisplacementsGCI'*BDisplacementsGCI)*elementLengthOnGP;
            end
            
            % Compute the rotations stiffness matrix
            if ~strcmp(alphaR,'undefined')
                KpRotI(EFTI,EFTI) = KpRotI(EFTI,EFTI) + ...
                    alphaR*(BRotationsI'*BRotationsI)*elementLengthOnGP;
            end
            
            % For patch J :
            % _____________
            
            % Compute the displacement stiffness matrix
            if ~strcmp(alphaD,'undefined')
                KpDispJ(EFTJ,EFTJ) = KpDispJ(EFTJ,EFTJ) + ...
                    alphaD*(BDisplacementsGCJ'*BDisplacementsGCJ)*elementLengthOnGP;
            end
            
            % Compute the rotations stiffness matrix
             if ~strcmp(alphaR,'undefined')
                KpRotJ(EFTJ,EFTJ) = KpRotJ(EFTJ,EFTJ) + ...
                    alphaR*(BRotationsJ'*BRotationsJ)*elementLengthOnGP;
             end
            
            %% 4iii.18. Compute the element coupling matrices at the GP
        
            % Compute the displacement coupling matrix
            if ~strcmp(alphaD,'undefined')
                CpDispI(EFTI,EFTJ) = CpDispI(EFTI,EFTJ) - ...
                    alphaD*(BDisplacementsGCI'*BDisplacementsGCJ)*elementLengthOnGP;
            end

            % Compute the rotation coupling matrix
            if ~strcmp(alphaR,'undefined')
                if haveSurfaceNormalsSameOrientation
                    CpRotI(EFTI,EFTJ) = CpRotI(EFTI,EFTJ) + ...
                        alphaR*(BRotationsI'*BRotationsJ)*elementLengthOnGP;
                else
                    CpRotI(EFTI,EFTJ) = CpRotI(EFTI,EFTJ) - ...
                        alphaR*(BRotationsI'*BRotationsJ)*elementLengthOnGP;
                end
            end
        end
    end
end

end
