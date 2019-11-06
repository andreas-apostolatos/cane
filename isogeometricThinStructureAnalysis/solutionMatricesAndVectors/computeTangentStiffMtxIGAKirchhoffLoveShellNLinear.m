function [tanStiffMtx,resVct,BSplinePatch,propCoupling,minElArea] = ...
    computeTangentStiffMtxIGAKirchhoffLoveShellNLinear...
    (constMtx,tanMtxLoad,dHat,dHatSaved,dHatDot,dDotSaved,BSplinePatch,...
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
% Returns the tangent stiffness matrix and internal force vector for the
% nonlinear analysis case to the isogeometrc Kirchhoff-Love shell.
%
%                 Input :
%              constMtx : The constant part of the tangent stiffness matrix
%                         and residual load vector (dummy variable for this 
%                         function)
%            tanMtxLoad : Tangent stiffness matrix resulting from the
%                         application of followe loads
%                  dHat : The displacement field of the previous iteration
%                         step (dummy variable for this function)
%             dHatSaved : The displacement field of the previous time step
%                         (dummy variable for this function)
%               dHatDot : The velocity field of the previous iteration step
%                         (dummy variable for this function)
%             dDotSaved : The velocity field of the previous time step 
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
%                         .DOFNumbering: Numbering of the DOFs sorted into
%                                        a 3D array
%                           .parameters: material parameters of the shell
%                                  .int: On the numerical integration
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
%               noPatch : Dummy array for this function
%            noTimeStep : Number of time step
%  noNonlinearIteration : Properties of the transient analysis module
%                         (dummy variable for this function)
%          noWeakDBCCnd : Number of weak Dirichlet boundary conditions
%                     t : The time instance
% propTransientAnalysis : Structure on the transient analysis properties :
%                               .timeDependence : 'true' or 'false'
%                                  .noTimeSteps : Number of time steps
%    isReferenceUpdated : Flag on whether the reference configuration is
%                         updated
%                   tab : Tabulation for outputting information onto the
%                         command window
%                outMsg : Allows outputting information onto the command
%                         window when chosen as 'outputEnabled'
%
%              Output :
%         tanStiffMtx : The master tangent stiffness matrix
%              resVct : The residual vector for the patch 
%                       resVct = FI + FExternal
%        BSplinePatch : The updated B-Spline patch cell
%        propCoupling : The updated with the stabilization terms coupling 
%                       properties array
%           minElArea : The minimum element area in the isogeometric mesh
%             
% Function layout :
%
% 0. Read input
%
% 1. Compute the material matrices
%
% 2. Choose an integration rule
%
% 3. loops over elements
% ->
%    3i. Create an element freedome table
%
%   3ii. Compute the determinant of the Jacobian to the transformation from the NURBS space (xi-eta) to the integration domain [-1,1]x[-1,1]
%
%  3iii. Initialize the element area
%
%   3iv. Loop over all Gauss points
%   ->
%        3iv.1. Compute the NURBS coordinates xi,eta of the Gauss Point coordinates in the bi-unit interval [-1, 1]
%
%        3iv.2. Compute the NURBS basis function and up to their second derivatives at the Gauss Point
%
%        3iv.3. Compute the covariant base vectors of the reference configuration and their first derivatives
%
%        3iv.4. Compute the surface normal of the reference configuration (third covariant base vector not normalized)
%
%        3iv.5. Compute the legth of G3Tilde (= area dA of the undeformed configuration)
%
%        3iv.6. Compute covariant base vectors of the current configuration and their first derivatives
%
%        3iv.7. Compute the surface normal of the reference configuration (third covariant base vector not normalized)
%
%        3iv.8. Compute element tangent stiffness matrix and residual internal load vector
%
%        3iv.9 Compute the element area on the Gauss Point and add the contribution
%
%        3iv.10. Add the contributions from the Gauss Point to the global matrix
%   <-
%    3v. Find the minimum element area in the isogemetric mesh
% <-
% 4. Assign the symmetric part to the stiffness matrix
%
% 5. Add existing tangent matrix from the application of follower loads
%
% 6. Compute the complete residual vector within the patch
%
% 7. Check output
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
int = BSplinePatch.int;

% Number of knots in xi-,eta-direction
mxi = length(Xi);
meta = length(Eta);
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));

% Check input
checkInputForBSplineSurface(p,mxi,nxi,q,meta,neta);

% Initialize minimum element area in the IGA mesh
tolerance = 1e-4;
if abs(CP(1,1,1)-CP(nxi,1,1)) >= tolerance
    minElArea = abs(CP(1,1,1)-CP(nxi,1,1));
else
    minElArea = CP(1,1,1)-CP(1,neta,1);
end

% Number of degrees of freedom for the whole structure
nDOFs = 3*nxi*neta;

% Number of degrees of freedom for the element
nDOFsEl = 3*(p+1)*(q+1);

% Initialize global stiffness matrix and internal force vector
tanStiffMtx = zeros(nDOFs,nDOFs);
FI = zeros(nDOFs,1);

%% 1. Compute the material matrices

% Compute the membrane material matrix
Dm = parameters.E*parameters.t/(1-parameters.nue^2)*...
    [1              parameters.nue 0
	 parameters.nue 1              0
     0              0              (1-parameters.nue)/2];
                                                 
% Compute the bending material matrix
Db = parameters.E*parameters.t^3/(12*(1-parameters.nue^2))*...
    [1              parameters.nue 0
     parameters.nue 1              0 
     0              0              (1-parameters.nue)/2];

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
        if (Xi(i+1)~=Xi(i) && Eta(j+1)~=Eta(j))
            %% 3i. Create an element freedome table
            
            % Initialize element freedome table
            EFT = zeros(1,nDOFsEl);
            
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
            
            %% 3iii. Initialize the element area
            elementArea = 0;
            
            %% 3iv. Loop over all Gauss points
            for cEta = 1:etaNGP
                for cXi = 1:xiNGP
                    %% 3iv.1. Compute the NURBS coordinates xi,eta of the Gauss Point coordinates in the bi-unit interval [-1, 1]
                    xi = ( Xi(i+1)+Xi(i) + xiGP(cXi)*(Xi(i+1)-Xi(i)) )/2;
                    eta = ( Eta(j+1)+Eta(j) + etaGP(cEta)*(Eta(j+1)-Eta(j)) )/2;
                    
                    %% 3iv.2. Compute the NURBS basis function and up to their second derivatives at the Gauss Point
                    nDrvBasis = 2;
                    dR = computeIGABasisFunctionsAndDerivativesForSurface...
                        (i,p,xi,Xi,j,q,eta,Eta,CP,isNURBS,nDrvBasis);
                    
                    %% 3iv.3. Compute the covariant base vectors of the reference configuration and their first derivatives
                    nDrvBaseVct = 1;
                    [dG1,dG2] = computeBaseVectorsAndDerivativesForBSplineSurface...
                        (i,p,j,q,CP,nDrvBaseVct,dR);
                 
                    %% 3iv.4. Compute the surface normal of the reference configuration (third covariant base vector not normalized)
                    G3Tilde = cross(dG1(:,1),dG2(:,1));
                    
                    %% 3iv.5. Compute the legth of G3Tilde (= area dA of the undeformed configuration)
                    dA = norm(G3Tilde);
                    
                    %% 3iv.6. Compute covariant base vectors of the current configuration and their first derivatives
                    nDrvBaseVctCur = 1;
                    [dg1,dg2] = computeBaseVectorsAndDerivativesForBSplineSurface...
                        (i,p,j,q,CPd,nDrvBaseVctCur,dR);
                    
                    %% 3iv.7. Compute the surface normal of the reference configuration (third covariant base vector not normalized)
                    g3Tilde = cross(dg1(:,1),dg2(:,1));
                    
                    %% 3iv.8. Compute element tangent stiffness matrix and residual internal load vector
                    [KTElOnGP,FIElOnGP] = computeIGAElTangentStiffMtxKirchhoffLoveShellNLinear...
                        (p,q,dR,[dG1(:,1) dG2(:,1)],[dG1(:,2) dG2(:,2) dG1(:,3)],...
                        G3Tilde,[dg1(:,1) dg2(:,1)],[dg1(:,2) dg2(:,2) dg1(:,3)],...
                        g3Tilde,Dm,Db);
                    
                    %% 3iv.9 Compute the element area on the Gauss Point and add the contribution
                    elementAreaOnGP = dA*detJxiu*xiGW(cXi)*etaGW(cEta);
                    elementArea = elementArea + elementAreaOnGP;
                    
                    %% 3iv.10. Add the contributions from the Gauss Point to the global matrix
                    tanStiffMtx(EFT,EFT) = tanStiffMtx(EFT,EFT) + KTElOnGP*elementAreaOnGP;
                    FI(EFT) = FI(EFT) + FIElOnGP*elementAreaOnGP;
                end
            end
            %% 3v. Find the minimum element area in the isogemetric mesh
            if elementArea<minElArea
                minElArea = elementArea;
            end
        end
    end
end

%% 4. Assign the symmetric part to the stiffness matrix
KTt = tanStiffMtx';
for i=1:nDOFs
    KTt(i,i) = 0;
end
tanStiffMtx = tanStiffMtx + KTt;

%% 5. Add existing tangent matrix from the application of follower loads
if ~ischar(tanMtxLoad)
    tanStiffMtx = tanStiffMtx + tanMtxLoad;
end

%% 6. Compute the complete residual vector within the patch
resVct = - (FI + loadFactor*BSplinePatch.FGamma);

%% 7. Check output
if isBSplinePatchCell
    BSplinePatch = {BSplinePatch};
end

end
