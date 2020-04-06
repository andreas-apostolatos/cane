function [K, F, BSplinePatches, propCoupling, minElAreaSize] = ...
    computeStiffMtxAndLoadVctDDMPenaltyIGAKirchoffLoveShell...
    (KConstantPenalty, tanMtxLoad, dHat, dHatSaved, dHatDot, dHatDotSaved, ...
    BSplinePatches, connections, propCoupling, loadFactor, noPatch, ...
    noTimeStep, iNLinearIter, noWeakDBCCnd, t, propStrDynamics, ...
    isReferenceUpdated, tab, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the coupled stiffness matrix of a multipatch system corresponding
% to the penalty formulation for the geometrically linear Kirchhoff-Love
% shell.
%
%              Input :
%   KConstantPenalty : The constant part of the coupled stiffness matrix
%                      using the Penalty method
%         tanMtxLoad : Dummy variable for this function
%               dHat : The complete displacement field of the multipatch 
%                      system
%          dHatSaved : Dummy variable for this function
%            dHatDot : Dummy variable for this function
%       dHatDotSaved : Dummy variable for this function
%     BSplinePatches : Array of B-Spline patches each of which containing 
%                      the following variables,
%                           .p,.q : The polynomial degrees
%                        .Xi,.Eta : The knot vectors
%                             .CP : The Control Poin coordinates and weights
%                        .isNURBS : Structure on whether the basis is a 
%                                   BSpline or a NURBS
%                            .int : Structure on the domain integration of 
%                                   the stiffness and body force vector
%                   .DOFNumbering : The numbering of the DOFs within the
%                                   patch itself
%        connections : Define the connection between the patches:
%                            .No : Number of connections
%                     .xiEtaCoup : [patchID1 patchID2 xi12 eta12 xi21 eta21
%                                   ...      ...    ...   ...  ...   ...]
%       propCoupling : Properties of the multipatch coupling
%                           .alphaD : penalty factor for the displacement
%                                     coupling
%                           .alphaR : penalty factor for the rotation
%                                     coupling
%                             .intC : On the integration of the coupling
%                                     interface
%         loadFactor : Dummy variable for this function
%            noPatch : Dummy variable for this function
%         noTimeStep : Dummy variable for this function
%       iNLinearIter : Dummy variable for this function
%       noWeakDBCCnd : Dummy variable for this function
%                  t : The time instance
%    propStrDynamics : Dummy variable for this function
% isReferenceUpdated : Dummy variable for this function
%                tab : The tabulation (dummy variable for this function)
%             outMsg : Enables output into the cell when its chosen as 
%                      'outputEnabled' (dummy variable for this function)
%
%             Output :
%                  K : The coupled with the penalty method stiffness matrix
%                  F : The externally applied load vector
%     BSplinePatches : Dummy output for this function
%       propCoupling : Dummy output for this function
%      minElAreaSize : Array containing the minimum element area sizes for
%                      each patch
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the patches of the multi-patch shell
% ->
%    1i. Compute the tangent stiffness matrix for the current patch
%
%   1ii. Assemble to the master tangent stiffness matrix and residual vector
% <-
%
% 2. Sum up the decoupled stiffness matrix with the penalty contributions
%
%% Function main body

%% 0. Read input

% Dummy variable
KConstantDummy = 'undefined';

% Compute the number of patches
numPatches = length(BSplinePatches);

% Compute the number of DOFs
numDOFs = length(dHat);

% Initialize the load vector
F = zeros(numDOFs,1);

% Initialize the stiffness matrix
K = zeros(numDOFs);

% Initialize the minimum element area size array
minElAreaSize = zeros(numPatches, 1);

%% 1. Loop over all the patches of the multi-patch shell
for iPatches = 1:numPatches
    %% 1i. Compute the tangent stiffness matrix for the current patch
    [KPatch, FPatch, ~, ~, minElAreaSize(iPatches,1)] = ...
        computeStiffMtxAndLoadVctIGAKirchhoffLoveShellLinear ...
        (KConstantDummy, tanMtxLoad, dHat, dHatSaved, dHatDot, ...
        dHatDotSaved, BSplinePatches{iPatches}, connections, ...
        propCoupling, loadFactor, noPatch, noTimeStep, iNLinearIter, ...
        noWeakDBCCnd, t, propStrDynamics, isReferenceUpdated, tab, ...
        outMsg);

    %% 1ii. Assemble to the master tangent stiffness matrix and residual vector
    K(BSplinePatches{iPatches}.EFTPatches,...
        BSplinePatches{iPatches}.EFTPatches) = KPatch;
    F(BSplinePatches{iPatches}.EFTPatches) = FPatch;
end

%% 2. Sum up the decoupled stiffness matrix with the penalty contributions
K = K + KConstantPenalty;

end