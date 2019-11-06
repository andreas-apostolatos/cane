function [errorInL2,errorInL2Gamma,exactNorm] = ...
    computeDomainAndInterfaceErrorInL2NormKLShellMultipatches...
    (BSplinePatches,dHat,connections,computeExRes,propProblem,outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the domain and interface jump error in the L2-norm corresponding 
% to the multipatch analysis of a Kirchhoff-Love shell as well as the 
% L2-norm of the exact displacement field in the domain.
%
%          Input :
% BSplinePatches : Structure containing all the information regarding
%                        the connections between the multipatches
%           dHat : The resulting displacement field distributed among the 
%                  multipatch geometry
%    connections : Define the connection between the patches:
%                            .No : Number of connections
%                     .xiEtaCoup : [patchID1 patchID2 xi12 eta12 xi21 eta21
%                                 ...      ...    ...   ...  ...   ...]
%   
%   computeExRes : Function handle to the computation of the exact
%                  resultants
%    propProblem : The characteristics of the problem for the computation
%                  of the analytical solution
%         outMsg : 'outputEnabled' enables output information into the
%                  command window while computations are taking place
%
%         Output :
%      errorInL2 : The L2-norm of the displacement error in the domain 
% errorInL2Gamma : The L2-norm of the interface jump in the displacement
%                  field
%      exactNorm : The L2-norm of the exact displacement field in the
%                  domain
%   
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the patches in the coupled system
% ->
%    1i. Recover the patch data
%
%   1ii. Get the discrete solution vector of the current patch
%
%  1iii. Select an integration scheme for the current patch
%
%   1iv. Distribute the discrete solution vector of the patch among the elements
%
%    1v. Loop over the elements of the patch
%    ->
%        1v.1. Compute the determinant of the Jacobian to the transformation from the NURBS space (xi-eta) to the integration domain [-1,1]x[-1,1] 
%
%        1v.2. Loop over the Gauss Points
%        ->
%              1v.2i. Compute the NURBS coordinates xi,eta of the Gauss Point coordinates in the bi-unit interval [-1, 1]
%
%             1v.2ii. Compute the NURBS basis function and up to their second derivatives at the Gauss Point
%
%            1v.2iii. Compute the Cartesian coordinates of the Gauss Point
%
%              1v.2v. Compute the covariant base vectors and their first derivatives
%
%             1v.2vi. Compute the surface normal (third covariant base vector not normalized)
%
%            1v.2vii. Compute the legth of G3Tilde (= area dA)
%
%           1v.2viii. Compute the exact resultants at the Gauss Point and transform them from the contravariant to the local Cartesian basis
%
%             1v.2ix. Compute the numerical displacement field and membrane strain energies of the error at the local Cartesian system at the Gauss Point
%
%              1v.2x. Compute the norm of the difference of the numerically computed and the analytical displacement fields in the L2-norm and add the contribution
%
%             1v.2xi. Compute the L2-norm of the analytical displacement field and add the contribution
%        <-
%
%    <-
%
% <-
%
% 2. loop over the connections between the patches to compute the interface jumps in the L2-norm
% ->
%    2i. Get the patch IDs
%
%   2ii. Get the patch parameters
%
%  2iii. Check if the patches have the same or opposite orientation
%
%   2iv. Select an integration scheme
%
%    2v. Get the discrete solution vector for each patch
%
%   2vi. Distribute the discrete solution vector for each patch among the elements
%
%  2vii. Get the running and the fixed parameters on the patch interface and the coupling region
%
% 2viii. Compute the merged knot vector from both patches over the interface
%
%  2ix. Loop over the elements on the coupling surface
%  ->
%       2ix.1. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
%
%       2ix.2. Loop over the Gauss points
%       ->
%              2ix.2i. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
%
%             2ix.2ii. Compute the NURBS basis functions
%
%            2ix.2iii. Compute the edge length on the Gauss Point
%
%             2ix.2iv. Compute the element discrete displacement vector for each patch
%
%              2ix.2v. Compute the displacement vector at the Gauss Point of the interface for both patches
%
%             2ix.2vi. Compute the difference of the displacements in the L2-norm on the Gauss Point and add the contribution
%       <-
%
%  <-
%
% <-
%
% 3. Take the square roots of the errors
%
% 4. Appendix 
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('____________________________________________________________________\n');
    fprintf('####################################################################\n');
    fprintf('Computation of the broken norm corresponding to the multipatch\n');
    fprintf('analysis of a Kirchhoff-Love shell has been initiated \n');
    fprintf('\n');
    fprintf('The employed example is a %s \n',propProblem.problem);
    if strcmp(propProblem.exactSolution,'analyticalSolution')
        fprintf('for which an analytical solution exists\n');
    elseif strcmp(propProblem.exactSolution,'overkillSolution')
        fprintf('which solution is compared against an overkill solution using \n one patch\n');
    else
        error('analytical solution to the problem not specified');
    end
    fprintf('____________________________________________________________________\n');
    fprintf('\n');
    
    % start measuring computational time
    tic;
end

%% 0. Read input

% Initialize output
errorInL2 = 0;
errorInL2Gamma = 0;
exactNorm = 0;

% Number of patches in the coupled system
noPatches = length(BSplinePatches);

% Create an element freedom tables for each patch
for i = 1:noPatches
    if i==1
        BSplinePatches{i}.EFTPatches = 1:3*BSplinePatches{i}.noCPs;
    else
        BSplinePatches{i}.EFTPatches = ...
            BSplinePatches{i-1}.EFTPatches(length(BSplinePatches{i-1}.EFTPatches)) + ...
            1:BSplinePatches{i-1}.EFTPatches(length(BSplinePatches{i-1}.EFTPatches)) + ...
            3*BSplinePatches{i}.noCPs;
    end
end

%% 1. Loop over all the patches in the coupled system
for counterPatches = 1:noPatches
    %% 1i. Recover the patch data
    p = BSplinePatches{counterPatches}.p;
    q = BSplinePatches{counterPatches}.q;
    Xi = BSplinePatches{counterPatches}.Xi;
    Eta = BSplinePatches{counterPatches}.Eta;
    CP = BSplinePatches{counterPatches}.CP;
    isNURBS = BSplinePatches{counterPatches}.isNURBS;
    int = BSplinePatches{counterPatches}.int;
    EFTPatch = BSplinePatches{counterPatches}.EFTPatches;
    mxi = length(Xi);
    meta = length(Eta);
    nxi = length(CP(:,1,1));
    
    %% 1ii. Get the discrete solution vector of the current patch
    dHatPatch = dHat(EFTPatch);
    
    %% 1iii. Select an integration scheme for the current patch
    if strcmp(int.type,'default')
        xiNGP = p + 1;
        etaNGP = q + 1;
    elseif strcmp(int.type,'user')
        xiNGP = int.xiNGP;
        etaNGP = int.etaNGP;
    end
    [xiGP,xiGW] = getGaussPointsAndWeightsOverUnitDomain(xiNGP);
    [etaGP,etaGW] = getGaussPointsAndWeightsOverUnitDomain(etaNGP);
    
    %% 1iv. Distribute the discrete solution vector of the patch among the elements
    dHatEl = zeros(mxi-p-1,meta-q-1,3*(p+1)*(q+1));
    for etaSpan = (q+1):(meta-q-1)
        for xiSpan = (p+1):(mxi-p-1)
            counterEl = 1; 
            for c = etaSpan-q-1:etaSpan-1 
                for b = xiSpan-p:xiSpan
                    dHatEl(xiSpan,etaSpan,counterEl) = dHatPatch(3*(c*nxi+b)-2);
                    dHatEl(xiSpan,etaSpan,counterEl + 1) = dHatPatch(3*(c*nxi+b)-1);
                    dHatEl(xiSpan,etaSpan,counterEl + 2) = dHatPatch(3*(c*nxi+b));

                    % Update counter
                    counterEl = counterEl + 3;
                end
            end
        end
    end
    if exist('dHatActual','var')
        clear dHatActual;
    end
      
    %% 1v. Loop over the elements of the patch
    for j = q+1:meta-q-1
        for i = p+1:mxi-p-1
            if Xi(i+1)~=Xi(i) && Eta(j+1)~=Eta(j)
                %% 1v.1. Compute the determinant of the Jacobian to the transformation from the NURBS space (xi-eta) to the integration domain [-1,1]x[-1,1] 
                %
                %         | xi_i+1 - xi_i                    |
                %         | -------------            0       |
                %         |        2                         |
                %  xi,u = |                                  |
                %         |                  eta_j+1 - eta_j |
                %         |        0         --------------- |
                %         |                          2       |
                detJxiu = (Xi(i+1)-Xi(i))*(Eta(j+1)-Eta(j))/4;
                
                %% 1v.2. Loop over the Gauss Points
                for cEta = 1:length(etaGW)
                    for cXi = 1:length(xiGW)
                        %% 1v.2i. Compute the NURBS coordinates xi,eta of the Gauss Point coordinates in the bi-unit interval [-1, 1]
                        xi = ( Xi(i+1)+Xi(i) + xiGP(cXi)*(Xi(i+1)-Xi(i)) )/2;
                        eta = ( Eta(j+1)+Eta(j) + etaGP(cEta)*(Eta(j+1)-Eta(j)) )/2;
                        
                        %% 1v.2ii. Compute the NURBS basis function and up to their second derivatives at the Gauss Point
                        nDrvBasis = 1;
                        dR = computeIGABasisFunctionsAndDerivativesForSurface...
                            (i,p,xi,Xi,j,q,eta,Eta,CP,isNURBS,nDrvBasis);
                        
                        %% 1v.2iii. Compute the Cartesian coordinates of the Gauss Point
                        P = computeCartesianCoordinatesOfAPointOnBSplineSurface...
                            (i,p,xi,Xi,j,q,eta,Eta,CP,dR(:,1));
                    
                        %% 1v.2v. Compute the covariant base vectors and their first derivatives
                        nDrvBaseVct = 0;
                        [A1,A2] = computeBaseVectorsAndDerivativesForBSplineSurface...
                            (i,p,j,q,CP,nDrvBaseVct,dR);

                        %% 1v.2vi. Compute the surface normal (third covariant base vector not normalized)
                        A3Tilde = cross(A1(:,1),A2(:,1));

                        %% 1v.2vii. Compute the legth of G3Tilde (= area dA)
                        dA = norm(A3Tilde);
                        
                        %% 1v.2viii. Compute the exact resultants at the Gauss Point and transform them from the contravariant to the local Cartesian basis
                        [d,epsilonRThetaContra,kappaRThetaContra] = ...
                            computeExRes(P(1,1),P(2,1),P(3,1),propProblem);
                        
                        %% 1v.2ix. Compute the numerical displacement field and membrane strain energies of the error at the local Cartesian system at the Gauss Point
                        dHatActual(:,1) = dHatEl(i,j,:);
                        dH = computePostprocDisplacementIGAKirchhoffLoveShell...
                            (p,q,dR(:,1),dHatActual);
                        
                        %% 1v.2x. Compute the norm of the difference of the numerically computed and the analytical displacement fields in the L2-norm and add the contribution
                        errorInL2 = errorInL2 + norm(d - dH)^2*dA*detJxiu*xiGW(cXi)*etaGW(cEta);
                        
                        %% 1v.2xi. Compute the L2-norm of the analytical displacement field and add the contribution
                        exactNorm = exactNorm + norm(d)^2*dA*detJxiu*xiGW(cXi)*etaGW(cEta);
                    end
                end
            end
        end
    end
end

%% 2. loop over the connections between the patches to compute the interface jumps in the L2-norm
for counterConnections = 1:connections.No
    %% 2i. Get the patch IDs
    
    % Patch 1 :
    % _________
    
    ID1 = connections.xiEtaCoup(counterConnections,1);
    
    % Patch 2 :
    % _________
    
    ID2 = connections.xiEtaCoup(counterConnections,2);
    
    %% 2ii. Get the patch parameters
    
    % Patch 1 :
    % _________
    
    p1 = BSplinePatches{ID1}.p;
    q1 = BSplinePatches{ID1}.q;
    Xi1 = BSplinePatches{ID1}.Xi;
    Eta1 = BSplinePatches{ID1}.Eta;
    CP1 = BSplinePatches{ID1}.CP;
    isNURBS1 = BSplinePatches{ID1}.isNURBS;
    EFTPatch1 = BSplinePatches{ID1}.EFTPatches;
    mxi1 = length(Xi1);
    meta1 = length(Eta1);
    nxi1 = length(CP1(:,1,1));
    neta1 = length(CP1(1,:,1));
    xicoup1 = connections.xiEtaCoup(counterConnections,3:4);
    etacoup1 = connections.xiEtaCoup(counterConnections,5:6);    
    
	% Patch 2 :
    % _________
    
    p2 = BSplinePatches{ID2}.p;
    q2 = BSplinePatches{ID2}.q;
    Xi2 = BSplinePatches{ID2}.Xi;
    Eta2 = BSplinePatches{ID2}.Eta;
    CP2 = BSplinePatches{ID2}.CP;
    isNURBS2 = BSplinePatches{ID2}.isNURBS;
    EFTPatch2 = BSplinePatches{ID2}.EFTPatches;
    mxi2 = length(Xi2);
    meta2 = length(Eta2);
    nxi2 = length(CP2(:,1,1));
    neta2 = length(CP2(1,:,1));
    xicoup2 = connections.xiEtaCoup(counterConnections,7:8);
    etacoup2 = connections.xiEtaCoup(counterConnections,9:10);
    
    %% 2iii. Check if the patches have the same or opposite orientation
    haveSameOrientation = findSubdomainInterfaceOrientation...
        (p1,Xi1,q1,Eta1,CP1,isNURBS1,xicoup1,etacoup1,p2,Xi2,q2,Eta2,CP2,isNURBS2,xicoup2,etacoup2);
    
    %% 2iv. Select an integration scheme
    p = max(p1,p2);
    q = max(q1,q2);
    pol = max(p,q);
    noGP = ceil((pol + 1)/2);
    [GP,GW] = getGaussPointsAndWeightsOverUnitDomain(noGP);
    
    %% 2v. Get the discrete solution vector for each patch
    
    % Patch 1 :
    % _________
    
    dHatPatch1 = dHat(EFTPatch1);
    
    % Patch 2 :
    % _________
    
    dHatPatch2 = dHat(EFTPatch2);
    
    %% 2vi. Distribute the discrete solution vector for each patch among the elements
    
    % Patch 1 :
    % _________
    
    dHatEl1 = zeros(mxi1-p1-1,meta1-q1-1,3*(p1+1)*(q1+1));
    for etaSpan = (q1+1):(meta1-q1-1)
        for xiSpan = (p1+1):(mxi1-p1-1)
            counterEl = 1; 
            for c = etaSpan-q1-1:etaSpan-1 
                for b = xiSpan-p1:xiSpan
                    dHatEl1(xiSpan,etaSpan,counterEl) = dHatPatch1(3*(c*nxi1+b)-2);
                    dHatEl1(xiSpan,etaSpan,counterEl + 1) = dHatPatch1(3*(c*nxi1+b)-1);
                    dHatEl1(xiSpan,etaSpan,counterEl + 2) = dHatPatch1(3*(c*nxi1+b));

                    % Update counter
                    counterEl = counterEl + 3;
                end
            end
        end
    end
    if exist('dHatActual1','var')
        clear dHatActual1;
    end
    
    % Patch 2 :
    % _________
    
    dHatEl2 = zeros(mxi2-p2-1,meta2-q2-1,3*(p2+1)*(q2+1));
    for etaSpan = (q2+1):(meta2-q2-1)
        for xiSpan = (p2+1):(mxi2-p2-1)
            counterEl = 1; 
            for c = etaSpan-q2-1:etaSpan-1 
                for b = xiSpan-p2:xiSpan
                    dHatEl2(xiSpan,etaSpan,counterEl) = dHatPatch2(3*(c*nxi2+b)-2);
                    dHatEl2(xiSpan,etaSpan,counterEl + 1) = dHatPatch2(3*(c*nxi2+b)-1);
                    dHatEl2(xiSpan,etaSpan,counterEl + 2) = dHatPatch2(3*(c*nxi2+b));

                    % Update counter
                    counterEl = counterEl + 3;
                end
            end
        end
    end
    if exist('dHatActual2','var')
        clear dHatActual2;
    end
    
    %% 2vii. Get the running and the fixed parameters on the patch interface and the coupling region

    % For patch 1 :
    % _____________

    if etacoup1(1) == etacoup1(2)
        % Coupled region in xi-direction
        couplingRegion1 = xicoup1;

        % Find the correct spans for the coupled region
        spanStart1 = findKnotSpan(couplingRegion1(1),Xi1,nxi1);
        spanEnd1 = findKnotSpan(couplingRegion1(2),Xi1,nxi1)+1;

        % Corresponding to the coupled region knot span
        couplingRegionOnKnotVector1 = Xi1(spanStart1:spanEnd1);

        % Fixed parameter on the parametric net
        eta1 = etacoup1(1);

        % Find the span where xiEta it lies in
        etaSpan1 = findKnotSpan(eta1,Eta1,neta1);

        % Flag on whether the coupling line is over xi
        isOnXi1 = true;
    else
        % Coupled region in eta-direction
        couplingRegion1 = etacoup1;

        % Find the correct spans for the coupled region
        spanStart1 = findKnotSpan(couplingRegion1(1),Eta1,neta1);   
        spanEnd1 = findKnotSpan(couplingRegion1(2),Eta1,neta1)+1;   

        % Corresponding to the coupled region knot span
        couplingRegionOnKnotVector1 = Eta1(spanStart1:spanEnd1);

        % Fixed parameter on the parametric net
        xi1 = xicoup1(1);

        % Find the span where uv it lies in
        xiSpan1 = findKnotSpan(xi1,Xi1,nxi1);

        % Flag on whether the coupling line is over eta
        isOnXi1 = false;
    end

    % For patch 2 :
    % _____________

    if etacoup2(1) == etacoup2(2)
        % Coupled region in xi-direction
        couplingRegion2 = xicoup2;

        % Find the correct spans for the coupled region
        spanStart2 = findKnotSpan(couplingRegion2(1),Xi2,nxi2);   
        spanEnd2 = findKnotSpan(couplingRegion2(2),Xi2,nxi2)+1; 

        % Corresponding to the coupled region knot span
        couplingRegionOnKnotVector2 = Xi2(spanStart2:spanEnd2);


        % Fixed parameter on the parametric net
        eta2 = etacoup2(1);

        % Find the span where xiEta it lies in
        etaSpan2 = findKnotSpan(eta2,Eta2,neta2);

        % Flag on whether the coupling line is over xi
        isOnXi2 = true;
    else
        % Coupled region in eta-direction
        couplingRegion2 = etacoup2;

        % Find the correct spans for the coupled region
        spanStart2 = findKnotSpan(couplingRegion2(1),Eta2,neta2);   
        spanEnd2 = findKnotSpan(couplingRegion2(2),Eta2,neta2)+1;

        % Corresponding to the coupled region knot span
        couplingRegionOnKnotVector2 = Eta2(spanStart2:spanEnd2);

        % Fixed parameter on the parametric net
        xi2 = xicoup2(1);

        % Find the span where uv it lies in
        xiSpan2 = findKnotSpan(xi2,Xi2,nxi2);

        % Flag on whether the coupling line is over eta
        isOnXi2 = false;
    end
    
    %% 2viii. Compute the merged knot vector from both patches over the interface
    knotVector = unique(mergesorted(couplingRegionOnKnotVector1,couplingRegionOnKnotVector2));
    
    %% 2ix. Loop over the elements on the coupling surface
    for i = 1:length(knotVector)-1
        if knotVector(i) ~= knotVector(i+1)
            %% 2ix.1. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
            detJxizeta = (knotVector(i+1)-knotVector(i))/2;
            
            %% 2ix.2. Loop over the Gauss points
            for j = 1:noGP
                %% 2ix.2i. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
                xiEta = ((1-GP(j))*knotVector(i)+(1+GP(j))*knotVector(i+1))/2;
                
                %% 2ix.2ii. Compute the NURBS basis functions
            
                % For patch 1 :
                % _____________

                if isOnXi1
                    xi1 = xiEta;
                    xiSpan1 = findKnotSpan(xi1,Xi1,nxi1);
                else
                    eta1 = xiEta;
                    etaSpan1 = findKnotSpan(eta1,Eta1,neta1);
                end
                dR1 = computeIGABasisFunctionsAndDerivativesForSurface...
                    (xiSpan1,p1,xi1,Xi1,etaSpan1,q1,eta1,Eta1,CP1,isNURBS1,1);

                % For patch 2 :
                % _____________

                if isOnXi2
                    xi2 = xiEta;
                    if ~haveSameOrientation
                        xi2 = Xi2(length(Xi2)) - xi2;
                    end
                    xiSpan2 = findKnotSpan(xi2,Xi2,nxi2);
                else
                    eta2 = xiEta;
                    if ~haveSameOrientation
                        eta2 = Eta2(length(Eta2)) - eta2;
                    end
                    etaSpan2 = findKnotSpan(eta2,Eta2,neta2);
                end
                R2 = computeIGABasisFunctionsAndDerivativesForSurface...
                    (xiSpan2,p2,xi2,Xi2,etaSpan2,q2,eta2,Eta2,CP2,isNURBS2,0);
                
                %% 2ix.2iii. Compute the edge length on the Gauss Point
                [A1,A2] = computeBaseVectorsAndDerivativesForBSplineSurface...
                    (xiSpan1,p1,etaSpan1,q1,CP1,0,dR1);
                if isOnXi1
                    dL = length(A1);
                else
                    dL = length(A2);
                end
                
                %% 2ix.2iv. Compute the element discrete displacement vector for each patch
                
                % For patch 1 :
                % _____________
                
                dHatActual1(:,1) = dHatEl1(xiSpan1,etaSpan1,:);
                
                % For patch 2 :
                % _____________
                
                dHatActual2(:,1) = dHatEl2(xiSpan2,etaSpan2,:);
                
                %% 2ix.2v. Compute the displacement vector at the Gauss Point of the interface for both patches
                
                % For patch 1 :
                % _____________
                
                dH1 = computePostprocDisplacementIGAKirchhoffLoveShell...
                    (p1,q1,dR1(:,1),dHatActual1);
                
                % For patch 2 :
                % _____________
                
                dH2 = computePostprocDisplacementIGAKirchhoffLoveShell...
                    (p2,q2,R2(:,1),dHatActual2);
                
                %% 2ix.2vi. Compute the difference of the displacements in the L2-norm on the Gauss Point and add the contribution
                errorInL2Gamma = errorInL2Gamma + norm(dH1 - dH2)^2*detJxizeta*dL*GW(j);
            end
        end
    end
end

%% 3. Take the square roots of the errors
errorInL2 = sqrt(errorInL2);
errorInL2Gamma = sqrt(errorInL2Gamma);
exactNorm = sqrt(exactNorm);
if strcmp(outMsg,'outputEnabled')
   fprintf('\t>> The L2-norm of the displacement error in the domain is %d\n',errorInL2);
   fprintf('\t>> The L2-norm of the displacement jump across the interface is %d\n',errorInL2Gamma);
end

%% 4. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;
    
    fprintf('Computation of the broken norm took %.2d seconds \n\n',computationalTime);
    fprintf('______________________Linear Analysis Ended_________________________\n');
    fprintf('####################################################################\n\n\n');
    fprintf('\n');
end

end
