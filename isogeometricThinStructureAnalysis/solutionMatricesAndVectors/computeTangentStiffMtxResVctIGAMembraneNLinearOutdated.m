function [tanStiffMtx,resVct,BSplinePatch,propCoupling,minElAreaSize] = ...
    computeTangentStiffMtxResVctIGAMembraneNLinearOutdated...
    (constMtx,tanMtxLoad,dHat,dHatSaved,dHatDot,dHatDotSaved,BSplinePatch,...
    connections,propCoupling,loadFactor,noPatch,noTimeStep,...
    noNonlinearIteration,noWeakDBCCnd,t,propTransientAnalysis,...
    isReferenceUpdated,tab,outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the tangent stiffness matrix and residual vector corresponding to
% the isogeometric membrane using the slow outdated algorithm.
%
%                 Input :
%              constMtx : The constant part of the tangent stiffness matrix
%                         from the application of weak boundary conditions
%                         using the penalty method
%            tanMtxLoad : Tangent stiffness matrix resulting from the
%                         application of followe loads
%                  dHat : The displacement field of the previous iteration
%                         step (dummy variable for this function)
%             dHatSaved : The displacement field of the previous time step
%                         (dummy variable for this function)
%               dHatDot : The velocity field of the previous iteration step
%                         (dummy variable for this function)
%          dHatDotSaved : The velocity field of the previous time step 
%                        (dummy variable for this function)
%          BSplinePatch : B-Spline patch with the following parameters:
%                                 .p,.q: Polynomial degrees
%                              .Xi,.Eta: knot vectors
%                                   .CP: Control Points coordinates and 
%                                        weights
%                              .isNURBS: Flag on whether the basis is a 
%                                        NURBS or a B-Spline
%                             .homDOFs : The global numbering of the
%                                        DOFs where homogeneous Dirichlet
%                                        boundary conditions are applied
%                           .inhomDOFs : The global numbering of the
%                                        DOFs where homogeneous Dirichlet
%                                        boundary conditions are applied
%                     .valuesInhomDOFs : Prescribed values to the DOFs 
%                                        where homogeneous Dirichlet
%                                        boundary conditions are applied
%                               FGamma : The boundary applied force vector
%                                        over the B-Spline patch
%                           bodyForces : Function handle to the computation 
%                                        of the body forces
%                        .DOFNumbering : Numbering of the DOFs sorted into
%                                        a 3D array
%                          .parameters : material parameters of the 
%                                        membrane
%                                 .int : On the numerical integration
%                                         .type : 'default' or 'user'
%                                        .xiNGP : No. of GPs along xi-
%                                                 direction for stiffness 
%                                                 entries
%                                       .etaNGP : No. of GPs along eta-
%                                                 direction for stiffness 
%                                                 entries
%                                 .xiNGPForLoad : No. of GPs along xi-
%                                                 direction for load 
%                                                 entries
%                                .etaNGPForLoad : No. of GPs along eta-
%                                                 direction for load 
%                                                 entries
%                                   .nGPForLoad : No. of GPs along boundary
%           connections : Dummy variable for this function
%          propCoupling : Dummy variable for this function
%            loadFactor : The load factor of the current time step
%               noPatch : Patch number in the multipatch array
%            noTimeStep : Number of time step
%  noNonlinearIteration : Number of the nonlinear iteration step
%          noWeakDBCCnd : Number of weak Dirichlet boundary conditions
%                     t : The time instance
% propTransientAnalysis : Structure on the transient analysis :
%                            .timeDependence : 'true' or 'false'
%                               .noTimeSteps : Number of time steps
%    isReferenceUpdated : Flag on whether the reference configuration is
%                         updated
%                   tab : Tabulation related to outputting information on
%                         the command window
%                outMsg : Enables outputting information onto the command
%                         window when chosen as 'outputEnabled'
%
%              Output :
%              tanMtx : The master tangent stiffness matrix
%              resVct : The global residual vector
%        BSplinePatch : The updated with the stabilization parameters
%                       BSpline patch array when weak enforcement of the
%                       Dirichlet boundary conditions is chosen
%        propCoupling : Updated Array of the coupling properties for
%                       multipatches coupled with the Nitsche method and
%                       automatic estimation of the stabilization is chosen
%       minElAreaSize : Dummy output for this function
%             
% Function layout :
%
% 0. Read input
%
% 1. Get the function handle for the computation of the tangent stiffness and the residual load vector corresponding to the application of weak Dirichlet boundary conditions
%
% 2. Choose an integration rule
%
% 3. loops over elements
% ->
%     3i. Create an element freedom table
%
%    3ii. Compute the determinant of the Jacobian to the transformation from the NURBS space (xi-eta) to the integration domain [-1,1]x[-1,1]
%
%   3iii. Loop over all Gauss points
%   ->
%         3iii.1. Compute the NURBS coordinates xi,eta of the Gauss Point coordinates in the bi-unit interval [-1, 1]
%
%         3iii.2. Compute the NURBS basis functions and their first derivatives at the Gauss Point
%
%         3iii.3. Compute the covariant base vectors of the reference configuration at the Gauss point
%
%         3iii.4. Compute the surface normal of the reference configuration at the Gauss point (third covariant base vector not normalized)
%
%         3iii.5. Compute the legth of G3Tilde at the Gauss point (= area dA of the undeformed configuration)
%
%         3iii.6. Compute covariant base vectors of the current configuration at the Gauss point
%
%         3iii.7. Compute element tangent stiffness matrix and residual internal load vector at the Gauss point
%
%         3iii.8. Compute the element area on the Gauss Point
%
%         3iii.9. Add the contributions from the Gauss Point to the global matrix
%   <-
% <-
%
% 4. Re-assemble the tangent stiffness matrix into a full matrix if only half of it is computed
%
% 5. Add existing tangent matrix from the application of follower loads
%
% 6. Update the residual load vector to account for the external loading
%
% 7. Updated the tangent stiffness matrix and the residual vector with the contributions resulting from the constant matrix
%
% 8. Compute the tangent matrix and residual vector contributions due to the application of weak boundary conditions
%
% 9. Compute the tangent stiffness matrix and residual load vector resulting from embedded into the patch cables
%
% 10. Check output
%
%% Function main body

%% 0. Read input

% Check the NURBS geometry input
isBSplinePatchCell = false;
if iscell(BSplinePatch)
    if length(BSplinePatch) > 1
        error('Multipatch NURBS surface is given as input to the computation of the stiffness matrix for a single patch NURBS surface');
    else
        BSplinePatch = BSplinePatch{1};
    end
    isBSplinePatchCell = true;
end

% Check the tangent matrix from the application of a follower load
if iscell(tanMtxLoad)
    if length(tanMtxLoad) > 1
        error('Multipatch tangent matrix due to follower loading is given as input to the computation of the stiffness matrix for a single patch NURBS surface');
    else
        tanMtxLoad = tanMtxLoad{1};
    end
end

% Get the patch parameters
p = BSplinePatch.p;
q = BSplinePatch.q;
Xi = BSplinePatch.Xi;
Eta = BSplinePatch.Eta;
CP = BSplinePatch.CP;
CPd = BSplinePatch.CPd;
isNURBS = BSplinePatch.isNURBS;
DOFNumbering = BSplinePatch.DOFNumbering;
parameters = BSplinePatch.parameters;
thickness = parameters.t; 
prestress = parameters.prestress;
int = BSplinePatch.int;

% Number of knots in xi-,eta-direction
mxi = length(Xi);
meta = length(Eta);
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));

% Assign dummy arrays
minElAreaSize = 'undefined';

% Check input
checkInputForBSplineSurface(p,mxi,nxi,q,meta,neta);

% Number of degrees of freedom for the whole structure
noDOFs = length(dHat);

% Number of degrees of freedom for the element
noDOFsEl = 3*(p+1)*(q+1);

% On the computation of the element tangent matrix and residual vector
% computeElTanMtxLoadVct = @computeIGAElTangentStiffMtxResVctMembraneNLinearOutdated;
computeElTanMtxLoadVct = @computeIGAElTangentStiffMtxResVctMembraneNLinearOutdatedKiendl;

% On the application of weak Dirichlet boundary conditions
isWeakDBC = false;
if ~isempty(BSplinePatch.weakDBC)
    if BSplinePatch.weakDBC.noCnd > 0
        isWeakDBC = true;
    end
end

% Initialize global stiffness matrix and internal force vector
tanStiffMtx = zeros(noDOFs,noDOFs);
if ~ischar(loadFactor)
    resVct = zeros(noDOFs,1);
else
    resVct = 'undefined';
end

% Assign the material matrix
Dm = parameters.E*parameters.t/(1-parameters.nue^2)*...
    [1              parameters.nue 0
	 parameters.nue 1              0
     0              0              (1-parameters.nue)/2];
 
%% 1. Get the function handle for the computation of the tangent stiffness and the residual load vector corresponding to the application of weak Dirichlet boundary conditions
if isWeakDBC
    if isfield(BSplinePatch.weakDBC,'method')
        if isfield(BSplinePatch.weakDBC,'imposedMotion')    
            if strcmp(BSplinePatch.weakDBC.method,'penalty')
                computeTangMtxResVct = @computeWeakDBCTangMtxResVctPenaltyIGAMembrane;
            elseif strcmp(BSplinePatch.weakDBC.method,'lagrangeMultipliers')
                computeTangMtxResVct = @computeWeakDBCTangMtxResVctLagrangeMultipliersIGAMembrane;
            end
        else
            if strcmp(BSplinePatch.weakDBC.method,'penalty') || ...
                    strcmp(BSplinePatch.weakDBC.method,'lagrangeMultipliers')
                computeTangMtxResVct = 'undefined';
            end
        end
        if strcmp(BSplinePatch.weakDBC.method,'nitsche')
            computeTangMtxResVct = @computeWeakDBCTangMtxResVctNitscheIGAMembrane;
        elseif ~strcmp(BSplinePatch.weakDBC.method,'penalty') && ...
               ~strcmp(BSplinePatch.weakDBC.method,'lagrangeMultipliers') && ...
               ~strcmp(BSplinePatch.weakDBC.method,'nitsche')
            error('Define a valid method for the application of weak Dirichlet boundary conditions in BSplinePatch.weakDBC.method')
        end
    else
        error('Define method for the application of weak Dirichlet boundary conditions in BSplinePatch.weakDBC');
    end
end
isWeakDBCTangent = false;
if isWeakDBC && isa(computeTangMtxResVct,'function_handle')
    isWeakDBCTangent = true;
end

%% 2. Choose an integration rule

% Select the integration scheme
if strcmp(int.type,'default')
    xiNGP = p + 1;
    etaNGP = q + 1;
elseif strcmp(int.type,'user')
    xiNGP = int.xiNGP;
    etaNGP = int.etaNGP;
end

% Issue the Gauss Point coordinates and weights
[xiGP,xiGW] = getGaussPointsAndWeightsOverUnitDomain(xiNGP);
[etaGP,etaGW] = getGaussPointsAndWeightsOverUnitDomain(etaNGP);

%% 3. loops over elements
for j = q+1:meta-q-1
    for i = p+1:mxi-p-1
        % check if element is greater than zero
        if Xi(i+1) ~= Xi(i) && Eta(j+1) ~= Eta(j)
            %% 3i. Create an element freedom table
            
            % Initialize element freedome table
            EFT = zeros(1,noDOFsEl);
            
            % initialize counter
            k = 1;
            
            % relation global-local dof
            for cpj = j-q:j
                for cpi = i-p:i
                    EFT(k) = DOFNumbering(cpi,cpj,1);
                    EFT(k+1) = DOFNumbering(cpi,cpj,2);
                    EFT(k+2) = DOFNumbering(cpi,cpj,3);
                    
                    % Update counter
                    k = k + 3;
                end
            end
            
            %% 3ii. Compute the determinant of the Jacobian to the transformation from the NURBS space (xi-eta) to the integration domain [-1,1]x[-1,1]
            %
            %         | xi_i+1 - xi_i                    |
            %         | -------------            0       |
            %         |        2                         |
            %  xi,u = |                                  |
            %         |                  eta_j+1 - eta_j |
            %         |        0         --------------- |
            %         |                          2       |
            detJxiu = (Xi(i+1)-Xi(i))*(Eta(j+1)-Eta(j))/4;
            
            %% 3iii. Loop over all Gauss points
            for cEta = 1:etaNGP
                for cXi = 1:xiNGP
                    %% 3iii.1. Compute the NURBS coordinates xi,eta of the Gauss Point coordinates in the bi-unit interval [-1, 1]
                    xi = ( Xi(i+1)+Xi(i) + xiGP(cXi)*(Xi(i+1)-Xi(i)) )/2;
                    eta = ( Eta(j+1)+Eta(j) + etaGP(cEta)*(Eta(j+1)-Eta(j)) )/2;
                    
                    %% 3iii.2. Compute the NURBS basis functions and their first derivatives at the Gauss Point
                    dR = computeIGABasisFunctionsAndDerivativesForSurface...
                        (i,p,xi,Xi,j,q,eta,Eta,CP,isNURBS,1);
                    
                    %% 3iii.3. Compute the covariant base vectors of the reference configuration at the Gauss point
                    [A1,A2] = computeBaseVectorsAndDerivativesForBSplineSurface...
                        (i,p,j,q,CP,0,dR);
                 
                    %% 3iii.4. Compute the surface normal of the reference configuration at the Gauss point (third covariant base vector not normalized)
                    A3Tilde = cross(A1(:,1),A2(:,1));
                    
                    %% 3iii.5. Compute the legth of G3Tilde at the Gauss point (= area dA of the undeformed configuration)
                    dA = norm(A3Tilde);
                    
                    %% 3iii.6. Compute covariant base vectors of the current configuration at the Gauss point
                    [a1,a2] = computeBaseVectorsAndDerivativesForBSplineSurface...
                        (i,p,j,q,CPd,0,dR);
                    
                    %% 3iii.7. Compute element tangent stiffness matrix and residual internal load vector at the Gauss point
                    [tanMtxElOnGP,resVctElOnGP] = computeElTanMtxLoadVct...
                        (i,p,xi,Xi,j,q,eta,Eta,CP,dR,[A1(:,1) A2(:,1)],...
                        [a1(:,1) a2(:,1)],thickness,Dm,prestress);
                    
                    %% 3iii.8. Compute the element area on the Gauss Point
                    elementAreaOnGP = dA*detJxiu*xiGW(cXi)*etaGW(cEta);
                    
                    %% 3iii.9. Add the contributions from the Gauss Point to the global matrix
                    tanStiffMtx(EFT,EFT) = tanStiffMtx(EFT,EFT) + tanMtxElOnGP*elementAreaOnGP;
                    if ~ischar(loadFactor)
                        resVct(EFT) = resVct(EFT) + resVctElOnGP*elementAreaOnGP;
                    end
                end
            end
        end
    end
end

%% 4. Re-assemble the tangent stiffness matrix into a full matrix if only half of it is computed
if strcmp(func2str(computeElTanMtxLoadVct),'computeIGAElTangentStiffMtxResVctMembraneNLinearOutdatedKiendl')
    tanMtxT = tanStiffMtx';
    for iCounter = 1:noDOFs
        tanStiffMtx(iCounter,iCounter) = 0;
    end
    tanStiffMtx = tanStiffMtx + tanMtxT;
end

%% 5. Add existing tangent matrix from the application of follower loads
if ~ischar(tanMtxLoad)
    tanStiffMtx = tanStiffMtx + tanMtxLoad;
end

%% 6. Update the residual load vector to account for the external loading
if ~ischar(loadFactor)
    resVct = resVct - loadFactor*BSplinePatch.FGamma;
end

%% 7. Updated the tangent stiffness matrix and the residual vector with the contributions resulting from the constant matrix
if ~ischar(constMtx)
    tanStiffMtx = tanStiffMtx + constMtx;
end
if ~ischar(constMtx) && ~ischar(loadFactor)
    resVct = resVct + constMtx*dHat;
end

%% 8. Compute the tangent stiffness matrix and residual load vector resulting from embedded into the patch cables
if isfield(BSplinePatch,'cables')
    if BSplinePatch.cables.No > 0
        [tanMtxCables,resVctCables] = computeTangentStiffMtxResVctCablesInThinStructureAnalysis...
            (BSplinePatch,noDOFs);
    end
else
    tanMtxCables = 'undefined';
    resVctCables = 'undefined';
end
if isfield(BSplinePatch,'cables')
    if BSplinePatch.cables.No > 0
        tanStiffMtx = tanStiffMtx + tanMtxCables;
    end
end
if isfield(BSplinePatch,'cables') && ~ischar(loadFactor)
    if BSplinePatch.cables.No > 0
        resVct = resVct + resVctCables;
    end
end

%% 9. Compute the tangent matrix and residual vector contributions due to the application of weak boundary conditions
if isWeakDBCTangent
    [tanMtxWeakDBC,resVctWeakDBC,BSplinePatch] = computeTangMtxResVct...
        (BSplinePatch,dHat,connections,noDOFs,propCoupling,tanStiffMtx,...
        noPatch,noTimeStep,noNonlinearIteration,noWeakDBCCnd,...
        thickness,t,propTransientAnalysis,tab,outMsg);
else
    tanMtxWeakDBC = 'undefined';
    resVctWeakDBC = 'undefined';
end
if ~ischar(tanMtxWeakDBC)
    tanStiffMtx = tanStiffMtx + tanMtxWeakDBC;
end
if ~ischar(resVctWeakDBC) && ~ischar(loadFactor)
    resVct = resVct + resVctWeakDBC;
end

%% 10. Check output
if isBSplinePatchCell
    BSplinePatch = {BSplinePatch};
end

end
