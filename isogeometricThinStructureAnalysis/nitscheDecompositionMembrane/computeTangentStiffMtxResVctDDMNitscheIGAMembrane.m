function [tanStiffMtx,resVct,BSplinePatches,propCoupling,minElAreaSize] = ...
    computeTangentStiffMtxResVctDDMNitscheIGAMembrane...
    (stabilMtx,tanMtxLoad,dHat,dHatSaved,dDotHat,dDotHatSaved,BSplinePatches,...
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
% Returns the tangent stiffness matrix and the residual vector
% corresponding to the nonlinear formulation of the Nitsche method for the
% decomposition of the isogeometric membrane.
%
%                 Input :
%             stabilMtx : The constant part of the tangent coupled system 
%                         matrix from the stabilization terms
%            tanMtxLoad : Tangent matrix resulting from the application of
%                         follower loading arranged in a cell array
%                  dHat : The complete control point displacement vector 
%                         from the previous iteration step
%             dHatSaved : The complete control point displacement vector 
%                         from the previous time step (dummy variable for 
%                         this function)
%               dDotHat : The complete control point velocity vector (dummy 
%                         variable for this function)
%          dDotHatSaved : The complete control point velocity vector from 
%                         the previous time step (dummy variable for this 
%                         function)
%        BSplinePatches : Its an array of structures {patch1,patch2,...} 
%                         each of the patch structures containing the 
%                         following information
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
%           connections : Define the connection between the patches:
%                            .No : Number of connections
%                     .xiEtaCoup : [patchID1 patchID2 xi12 eta12 xi21 eta21
%                                   ...      ...    ...   ...  ...   ...]
%          propCoupling : Properties of the multipatch coupling
%                           .alphaD : penalty factor for the displacement
%                                     coupling
%                           .alphaR : penalty factor for the rotation
%                                     coupling
%                             .intC : On the integration of the coupling
%                                     interface 
%            loadFactor : The load factor correspoding to the current load 
%                         step
%               noPatch : Dummy array for this function
%            noTimeStep : Number of time step
%  noNonlinearIteration : Number of the nonlinear iteration step
%          noWeakDBCCnd : Number of weak Dirichlet boundary conditions
%                     t : The time instance
% propTransientAnalysis : Structure on the transient analysis :
%                            .timeDependence : 'true' or 'false'
%                               .noTimeSteps : Number of time steps
%    isReferenceUpdated : Flag on whether the reference geometry is updated 
%                   tab : Tabulation for outputting information onto the
%                         command window
%                outMsg : Allowing outputting information onto the command
%                         window when chosen as 'outputEnabled'
%
%                Output :
%              KTangent : The tangent stiffness matrix corresponding to the
%                         multipatch coupling within the penalty method
%                ResVct : The residual vector corresponding to the
%                         multipatch coupling within the penalty method
%        BSplinePatches : The updated cell array of B-Spline patches
%          propCoupling : The updated with the stabilization terms
%                         coupling properties array
%         minElAreaSize : Array noPatches x 1 containing the minimum
%                         element area size for each patch
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the patches of the multipatch membrane
% ->
%    1i. Assemble to the external load vector for the current patch
%
%   1ii. Compute the tangent stiffness matrix for the current patch
%
%  1iii. Assemble to the master tangent stiffness matrix and residual vector
% <-
%
% 2. Loop over all the connections in the multipatches and assemble the global penalty matrix
% ->
%    2i. Get the patch IDs
%
%   2ii. Get the coupling boundaries for both patches
%
%  2iii. Determine the interface orientation
%
%   2iv. Compute the penalty contributions to the coupled system
%
%    2v. Assemble the Nitsche matrices of each patch to the master Nitsche stiffness/coupling matrix
%
%   2vi. Assemble the Nitsche residual vectors of each patch to the master Nitsche residual vectors
% <-
%
% 3. Compute the tangent stiffness matrix of the coupled system
%
% 4. Compute the residual vector of the coupled system
%
%% Function main body

%% 0. Read input

% Initialize dummy arrays
dHatDot = 'undefined';
dDotSaved = 'undefined';

% Initialize the structure array of the tangent stiffness matrix from each 
% patch needed when estimation of the stabilization is enabled
tanStiffMtxPatch = struct([]);

% Get the number of B-Spline patches
noPatches = length(BSplinePatches);

% Compute the number of DOFs
noDOFs = length(dHat);

% Initialize the master Nitsche matrix of the multipatch system
tanStiffMtxNitsche = zeros(noDOFs,noDOFs);

% Initialize the externally applied load vector
if ~ischar(loadFactor)
    FExternal = zeros(noDOFs,1);
else
    FExternal = 'undefined';
end

% Initialize the master Nitsche residual vector of the multipatch system
if ~ischar(loadFactor)
    resVctNitsche = zeros(noDOFs,1);
else
    resVctNitsche = 'undefined';
end

% Initialize master tangent stiffness matrix of the decoupled multipatch 
% system
tanStiffMtxDecoupled = zeros(noDOFs,noDOFs);

% Initialize the master residual vector of the multi-patch system
if ~ischar(loadFactor)
    resVct = zeros(noDOFs,1);
else
    resVct = 'undefined';
end

% Initialize the minimum element area size array
minElAreaSize = struct([]);

%% 1. Loop over all the patches of the multipatch membrane
for iPatches = 1:noPatches
    %% 1i. Assemble to the external load vector for the current patch
    if ~ischar(loadFactor)
        FExternal(BSplinePatches{iPatches}.EFTPatches) = ...
                loadFactor*BSplinePatches{iPatches}.FGamma;
    end

    %% 1ii. Compute the tangent stiffness matrix for the current patch
    if BSplinePatches{iPatches}.noElmnts == 1
        computeTangentStiffMtxResVct = @computeTangentStiffMtxResVctIGAMembraneNLinearOutdated;
    else
        computeTangentStiffMtxResVct = @computeTangentStiffMtxResVctIGAMembraneNLinear;
    end
    [tanStiffMtxPatch{iPatches},resVctPatch,BSplinePatches{iPatches},...
        propCoupling,minElAreaSize{iPatches}] = ...
        computeTangentStiffMtxResVct...
        (BSplinePatches{iPatches}.KConstant,tanMtxLoad{iPatches},...
        dHat(BSplinePatches{iPatches}.EFTPatches),dHatSaved,dHatDot,...
        dDotSaved,BSplinePatches{iPatches},connections,propCoupling,...
        loadFactor,iPatches,noTimeStep,noNonlinearIteration,...
        noWeakDBCCnd,t,propTransientAnalysis,isReferenceUpdated,tab,outMsg);
        
    %% 1iii. Assemble to the master tangent stiffness matrix and residual vector
    tanStiffMtxDecoupled(BSplinePatches{iPatches}.EFTPatches,...
        BSplinePatches{iPatches}.EFTPatches) = tanStiffMtxPatch{iPatches};
    if ~ischar(loadFactor)
        resVct(BSplinePatches{iPatches}.EFTPatches) = resVctPatch;
    end
end

%% 2. Loop over all the connections in the multipatches and assemble the global penalty matrix
for iConnections = 1:connections.No
    %% 2i. Get the patch IDs
    
    % Patch 1 :
    % _________
    
    ID1 = connections.xiEtaCoup(iConnections,1);
    
    % Patch 2 :
    % _________
    
    ID2 = connections.xiEtaCoup(iConnections,2);
    
    %% 2ii. Get the coupling boundaries for both patches
    
    % Patch 1 :
    % _________
    
    BSplinePatches{ID1}.xicoup = connections.xiEtaCoup(iConnections,3:4);
    BSplinePatches{ID1}.etacoup = connections.xiEtaCoup(iConnections,5:6);
    
	% Patch 2 :
    % _________
    
    BSplinePatches{ID2}.xicoup = connections.xiEtaCoup(iConnections,7:8);
    BSplinePatches{ID2}.etacoup = connections.xiEtaCoup(iConnections,9:10);
    
    %% 2iii. Determine the interface orientation
    haveSameDirection = findSubdomainInterfaceOrientation...
        (BSplinePatches{ID1}.p,BSplinePatches{ID1}.Xi,BSplinePatches{ID1}.q,BSplinePatches{ID1}.Eta,BSplinePatches{ID1}.CP,BSplinePatches{ID1}.isNURBS,BSplinePatches{ID1}.xicoup,BSplinePatches{ID1}.etacoup,...
        BSplinePatches{ID2}.p,BSplinePatches{ID2}.Xi,BSplinePatches{ID2}.q,BSplinePatches{ID2}.Eta,BSplinePatches{ID2}.CP,BSplinePatches{ID2}.isNURBS,BSplinePatches{ID2}.xicoup,BSplinePatches{ID2}.etacoup);
    
    %% 2iv. Compute the penalty contributions to the coupled system
    [KNitsche1,CNitsche1,KNitsche2,resVctNitsche1,resVctNitsche2,propCoupling] = ....
        computeDDMTangentNitscheMtcesIGAMembrane...
        (BSplinePatches{ID1},BSplinePatches{ID2},dHat(BSplinePatches{ID1}.EFTPatches),...
        dHat(BSplinePatches{ID2}.EFTPatches),haveSameDirection,connections,...
        propCoupling,tanStiffMtxPatch{ID1},tanStiffMtxPatch{ID2},...
        iConnections,noTimeStep,noNonlinearIteration,propTransientAnalysis,...
        tab,outMsg);
    
    % We can use the transpose of C1PenaltyDisplacements and
    % C1PenaltyRotations see below the assembly to the global coupling 
    % matrix for both patches
    
    %% 2v. Assemble the Nitsche matrices of each patch to the master Nitsche stiffness/coupling matrix
    
    % Assemble the stiffness matrix contribution to patch 1
    tanStiffMtxNitsche(BSplinePatches{ID1}.EFTPatches,BSplinePatches{ID1}.EFTPatches) = ...
        tanStiffMtxNitsche(BSplinePatches{ID1}.EFTPatches,BSplinePatches{ID1}.EFTPatches) + ...
        KNitsche1;
    
    % Assemble the stiffness matrix contribution to patch 2
    tanStiffMtxNitsche(BSplinePatches{ID2}.EFTPatches,BSplinePatches{ID2}.EFTPatches) = ...
        tanStiffMtxNitsche(BSplinePatches{ID2}.EFTPatches,BSplinePatches{ID2}.EFTPatches) + ...
        KNitsche2;
    
    % Assemble to the coupling matrix contribution of the variational
    % problem in patch 1
    tanStiffMtxNitsche(BSplinePatches{ID1}.EFTPatches,BSplinePatches{ID2}.EFTPatches) = ...
        tanStiffMtxNitsche(BSplinePatches{ID1}.EFTPatches,BSplinePatches{ID2}.EFTPatches) + ...
        CNitsche1;
    
    % Assemble to the coupling matrix contribution of the variational
    % problem in patch 2
    tanStiffMtxNitsche(BSplinePatches{ID2}.EFTPatches,BSplinePatches{ID1}.EFTPatches) = ...
        tanStiffMtxNitsche(BSplinePatches{ID2}.EFTPatches,BSplinePatches{ID1}.EFTPatches) + ...
        CNitsche1';
    
    %% 2vi. Assemble the Nitsche residual vectors of each patch to the master Nitsche residual vectors
    
    % Assemble the residual vector contribution to patch 1
    if ~ischar(loadFactor)
        resVctNitsche(BSplinePatches{ID1}.EFTPatches,1) = ...
            resVctNitsche(BSplinePatches{ID1}.EFTPatches,1) + ...
            resVctNitsche1;
    end
    
    % Assemble the residual vector contribution to patch 1
    if ~ischar(loadFactor)
        resVctNitsche(BSplinePatches{ID2}.EFTPatches,1) = ...
            resVctNitsche(BSplinePatches{ID2}.EFTPatches,1) + ...
            resVctNitsche2;
    end
end

%% 3. Compute the tangent stiffness matrix of the coupled system
if ~ischar(stabilMtx)
    tanStiffMtx = tanStiffMtxDecoupled + tanStiffMtxNitsche + stabilMtx;
else
    tanStiffMtx = tanStiffMtxDecoupled + tanStiffMtxNitsche;
end

%% 4. Compute the residual vector of the coupled system
if ~ischar(loadFactor)
    if ~ischar(stabilMtx)
        resVct = resVct + stabilMtx*dHat + resVctNitsche;
    else
        resVct = resVct + resVctNitsche;
    end
end

end
