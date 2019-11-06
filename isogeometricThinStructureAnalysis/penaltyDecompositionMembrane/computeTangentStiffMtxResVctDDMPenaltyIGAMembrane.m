function [tanStiffMtx,resVct,BSplinePatches,propCoupling,minElAreaSize] = ...
    computeTangentStiffMtxResVctDDMPenaltyIGAMembrane...
    (stiffMtxLM,tanMtxLoad,dHat,dHatSaved,dDotHat,dDotHatSaved,BSplinePatches,...
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
% corresponding to the nonlinear formulation of the penalty method for the
% decomposition of the isogeometric membrane.
%
%                 Input :
%            stiffMtxLM : The constant part of the tangent system matrix 
%                         for the coupled system using the augmented
%                         Lagrange multipliers method
%            tanMtxLoad : The tangent matrix resulting from the
%                         application of follower loading arranged in a
%                         cell array
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
%           connections : Dummy array for this function
%          propCoupling : Dummy array for this function
%            loadFactor : The load factor correspoding to the current load 
%                         step
%               noPatch : Dummy array for this function
%            noTimeStep : Number of time step
%  noNonlinearIteration : Number of the nonlinear iteration step
%          noWeakDBCCnd : Number of weak Dirichlet boundary conditions
%                     t : The time instance
% propTransientAnalysis : Structure for the transient analysis properties
%                               .timeDependence : 'true' or 'false'
%                                  .noTimeSteps : Number of time steps
%    isReferenceUpdated : Flag on whether the reference configuration is
%                         updated
%                   tab : Tabulation for outputting information on the
%                         command window
%                outMsg : Allows outputting information onto the command
%                         window when chosen as 'outputEnabled'
%
%                Output :
%           tanStiffMtx : The tangent stiffness matrix corresponding to the
%                         multipatch coupling within the penalty method
%                resVct : The residual vector corresponding to the
%                         multipatch coupling within the penalty method
%        BSplinePatches : The updated cell array of the B-Spline patches
%          propCoupling : The updated with the stabilization terms coupling 
%                         properties array
%         minElAreaSize : Array noPatches x 1 containing the minimum
%                         element area size for each patch
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the patches of the multipatch mebrane
% ->
%    1i. Assemble to the external load vector for the current patch
%
%   1ii. Compute the tangent stiffness matrix for the current patch
%
%  1iii. Assemble to the master tangent stiffness matrix and residual vector
% <-
% 
% 2. Compute the tangent stiffness matrix of the coupled system
%
% 3. Compute the residual vector of the coupled system
% 
%% Function main body

%% 0. Read input

% Initialize dummy arrays
dHatDot = 'undefined';
dDotSaved = 'undefined';

% Get the number of B-Spline patches
noPatches = length(BSplinePatches);

% Compute the number of DOFs
noDOFs = length(dHat);

% Initialize the externally applied load vector
FExternal = zeros(noDOFs,1);

% Initialize master tangent stiffness matrix of the multipatch system
tanStiffMtxDecoupled = zeros(noDOFs,noDOFs);

% Initialize the master residual vector of the multi-patch system
if ~ischar(loadFactor)
    resVctDecoupled = zeros(noDOFs,1);
end

% Initialize the minimum element area size array
minElAreaSize = struct([]);

%% 1. Loop over all the patches of the multi-patch mebrane
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
    [tanStiffMtxPatch,resVctPatch,BSplinePatches{iPatches},minElAreaSize{iPatches}] = ...
        computeTangentStiffMtxResVct...
        (BSplinePatches{iPatches}.KConstant,tanMtxLoad{iPatches},...
        dHat(BSplinePatches{iPatches}.EFTPatches),dHatSaved,dHatDot,...
        dDotSaved,BSplinePatches{iPatches},connections,propCoupling,...
        loadFactor,iPatches,noTimeStep,noNonlinearIteration,...
        noWeakDBCCnd,t,propTransientAnalysis,isReferenceUpdated,tab,outMsg);

    %% 1iii. Assemble to the master tangent stiffness matrix and residual vector
    tanStiffMtxDecoupled(BSplinePatches{iPatches}.EFTPatches,BSplinePatches{iPatches}.EFTPatches) = ...
        tanStiffMtxPatch;
    if ~ischar(loadFactor)
        resVctDecoupled(BSplinePatches{iPatches}.EFTPatches) = resVctPatch;
    end
end

%% 2. Compute the tangent stiffness matrix of the coupled system
if strcmp(propCoupling.method,'penalty') || strcmp(propCoupling.method,'lagrangeMultipliers')
    tanStiffMtx = tanStiffMtxDecoupled + stiffMtxLM;
elseif strcmp(propCoupling.method,'mortar')
    tanStiffMtx = tanStiffMtxDecoupled;
else
    error('This function call requires the coupling method to be either Penalty, Lagrange Multipliers, or Mortar but not %s',propCoupling.method);
end

%% 3. Compute the residual vector of the coupled system
if ~ischar(loadFactor)
    if strcmp(propCoupling.method,'penalty') || strcmp(propCoupling.method,'lagrangeMultipliers')
        resVct = resVctDecoupled + stiffMtxLM*dHat;
    elseif strcmp(propCoupling.method,'mortar')
        resVct = resVctDecoupled;
    end
else
    resVct = 'undefined';
end

end
