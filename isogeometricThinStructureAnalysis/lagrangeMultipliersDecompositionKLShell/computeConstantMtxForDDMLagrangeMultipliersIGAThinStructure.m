function KLagrangeMultipliers = computeConstantMtxForDDMLagrangeMultipliersIGAThinStructure...
    (BSplinePatches,connections,noDOFs,propCoupling)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the part of the tangent stiffness matrix which remains constant
% throughout the nonlinear iterations for the decomposition of the
% isogeometric membrane/Kirchhoff-Love shell using the Lagrange Multipliers 
% method.
%
%                Input :
%       BSplinePatches : Structure containing all the information regarding 
%                        the connections between the multipatches
%          connections : Define the connection between the patches:
%                            .No : Number of connections
%                     .xiEtaCoup : [patchID1 patchID2 xi12 eta12 xi21 eta21
%                                 ...      ...    ...   ...  ...   ...]
%                        .lambda : The NURBS discretization of the 
%                                  interface traction force field for each 
%                                  defined interface
%                            .mu : The NURBS discretization of the 
%                                  interface traction moment field for each 
%                                  defined interface
%               noDOFs : The complete number of DOFs
%         propCoupling : Properties of the multipatch coupling
%                           .alphaD : Vector containing the penalty factor
%                                     for the displacements and for each 
%                                     patch pair
%                           .alphaR : Vector containing the penalty factor
%                                     for the rotations and for each 
%                                     patch pair
%                             .intC : On the integration of the coupling
%                                     interface 
%
%               Output :
% KLagrangeMultipliers : The complete Lagrange Multipliers/penalty 
%                        contribution to the coupled system:
%
%                                   | Kp1      Cp1      Lambda1 Mu1  |
%           KLagrangeMultipliers =  | Kp2      Cp2      Lambda2 Mu2  |
%                                   | Lambda1' Lambda2' null    null |
%                                   | Mu1'     Mu2      null    null |
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the connections in the multipatches and assemble the global Lagrange Multipliers matrix
% ->
%    1i. Get the patch IDs
%
%   1ii. Get the penalty parameters for the current patch pair
%
%  1iii. Get the coupling boundaries for both patches
%
%   1iv. Get the Lagrange Multipliers fields for the current interface
%
%    1v. Determine the interface orientation
%
%   1vi. Compute the penalty contributions to the coupled system 
%
%  1vii. Compute the Lagrange Multipliers matrices
%
% 1viii. Assemble the penalty matrices to the master penalty stiffness/coupling matrix
%
%   1ix. Assemble the Lagrange Multipliers matrices to master saddle point matrix corrssponding to the variational forms of each patch
%
%   1xi. Assemble the Lagrange Multipliers matrices to master saddle point matrix corrssponding to the variational forms each Lagrange Multipliers field
% <-
%
%% Function main body

%% 0. Read input

% Initialize the master penalty matrix of the multi-patch system
% if strcmp(propCoupling.method,'lagrangeMultipliers')
%     KLagrangeMultipliers = zeros(noDOFs,noDOFs);
% elseif strcmp(propCoupling.method,'mortar')
%     KLagrangeMultipliers = struct([]);
% else
%     error('For this function the coupling method can be either lagrangeMultipliers or mortar but not %s',propCoupling.method);
% end
KLagrangeMultipliers = zeros(noDOFs,noDOFs);

% Assign the master and slave relation
isMaster = 1;
isSlave = 0;

% Check if Lagrange Multipliers field for the rotational coupling is
% enforced
isLagrangeMultipliersRotationsEnabled = false;
if isfield(connections,'mu')
    isLagrangeMultipliersRotationsEnabled = true;
    if isempty(connections.mu)
        error('Lagrange multipliers discretization for the rotational coupling exists but its empty');
    end
end

% Get the interface quadrature
intC = propCoupling.intC;

%% 1. Loop over all the connections in the multipatches and assemble the global Lagrange Multipliers matrix
for iConnections = 1:connections.No
    %% 1i. Get the patch IDs
    
    % Patch I :
    % _________
    
    idI = connections.xiEtaCoup(iConnections,1);
    
    % Patch J :
    % _________
    
    idJ = connections.xiEtaCoup(iConnections,2);
    
    %% 1ii. Get the penalty parameters for the current patch pair
    if isfield(propCoupling,'alphaD')
        alphaDIJ = propCoupling.alphaD(iConnections,1);
    else
        alphaDIJ = 0;
    end
    if isfield(propCoupling,'alphaR')
        alphaRIJ = propCoupling.alphaR(iConnections,1);
    else
        alphaRIJ = 0;
    end
    
    %% 1iii. Get the coupling boundaries for both patches
    
    % Patch I :
    % _________
    
    BSplinePatches{idI}.xicoup = connections.xiEtaCoup(iConnections,3:4);
    BSplinePatches{idI}.etacoup = connections.xiEtaCoup(iConnections,5:6);
    
	% Patch J :
    % _________
    
    BSplinePatches{idJ}.xicoup = connections.xiEtaCoup(iConnections,7:8);
    BSplinePatches{idJ}.etacoup = connections.xiEtaCoup(iConnections,9:10);
    
    %% 1iv. Get the Lagrange Multipliers fields for the current interface
    lambda = connections.lambda{iConnections};
    if isLagrangeMultipliersRotationsEnabled
        mu = connections.mu{iConnections};
    else
        mu = 'undefined';
    end
    
    %% 1v. Determine the interface orientation
    haveSameDirection = findSubdomainInterfaceOrientation...
        (BSplinePatches{idI}.p,BSplinePatches{idI}.Xi,BSplinePatches{idI}.q,BSplinePatches{idI}.Eta,BSplinePatches{idI}.CP,BSplinePatches{idI}.isNURBS,BSplinePatches{idI}.xicoup,BSplinePatches{idI}.etacoup,...
        BSplinePatches{idJ}.p,BSplinePatches{idJ}.Xi,BSplinePatches{idJ}.q,BSplinePatches{idJ}.Eta,BSplinePatches{idJ}.CP,BSplinePatches{idJ}.isNURBS,BSplinePatches{idJ}.xicoup,BSplinePatches{idJ}.etacoup);
    
    %% 1vi. Compute the penalty contributions to the coupled system 
    [KPenaltyDisplacementsI,KPenaltyRotationsI,CPenaltyDisplacementsI,...
        CPenaltyRotationsI,KPenaltyDisplacementsJ,KPenaltyRotationsJ] = ....
        computeDDMPenaltyMtcesIGAThinStructure...
        (BSplinePatches{idI},BSplinePatches{idJ},alphaDIJ,alphaRIJ,haveSameDirection,intC);
    
    % It can be used the transpose of C1PenaltyDisplacements and
    % C1PenaltyRotations see below the assembly to the global coupling 
    % matrix for both patches
    
    %% 1vii. Compute the Lagrange Multipliers matrices
    
    % Patch I :
    % _________
    
    [LambdaI,MuI] = computeDDMLagrangeMultipliersMtces4IGAKLShellLinear(BSplinePatches{idI},...
        lambda,mu,isMaster,haveSameDirection,propCoupling);
    
    % Patch J :
    % _________
    
    [LambdaJ,MuJ] = computeDDMLagrangeMultipliersMtces4IGAKLShellLinear(BSplinePatches{idJ},...
        lambda,mu,isSlave,haveSameDirection,propCoupling);
    
    %% 1viii. Assemble the penalty matrices to the master penalty stiffness/coupling matrix
    
    % Assemble the stiffness matrix contribution to patch I
    KLagrangeMultipliers(BSplinePatches{idI}.EFTPatches,BSplinePatches{idI}.EFTPatches) = ...
        KLagrangeMultipliers(BSplinePatches{idI}.EFTPatches,BSplinePatches{idI}.EFTPatches) + ...
        KPenaltyDisplacementsI + KPenaltyRotationsI;
    
    % Assemble the stiffness matrix contribution to patch J
    KLagrangeMultipliers(BSplinePatches{idJ}.EFTPatches,BSplinePatches{idJ}.EFTPatches) = ...
        KLagrangeMultipliers(BSplinePatches{idJ}.EFTPatches,BSplinePatches{idJ}.EFTPatches) + ...
        KPenaltyDisplacementsJ + KPenaltyRotationsJ;
    
    % Assemble to the coupling matrix contribution of the variational
    % problem in patch I
    KLagrangeMultipliers(BSplinePatches{idI}.EFTPatches,BSplinePatches{idJ}.EFTPatches) = ...
        KLagrangeMultipliers(BSplinePatches{idI}.EFTPatches,BSplinePatches{idJ}.EFTPatches) + ...
        CPenaltyDisplacementsI + CPenaltyRotationsI;
    
    % Assemble to the coupling matrix contribution of the variational
    % problem in patch J
    KLagrangeMultipliers(BSplinePatches{idJ}.EFTPatches,BSplinePatches{idI}.EFTPatches) = ...
        KLagrangeMultipliers(BSplinePatches{idJ}.EFTPatches,BSplinePatches{idI}.EFTPatches) + ...
        CPenaltyDisplacementsI' + CPenaltyRotationsI';
    
    %% 1ix. Assemble the Lagrange Multipliers matrices to master saddle point matrix corrssponding to the variational forms of each patch
    
    % Assemble the Lagrange multipliers field for the traction forces  of
    % the variational problem in patch I
    KLagrangeMultipliers(BSplinePatches{idI}.EFTPatches,lambda.EFTLagrangeMultipliers) = ...
        KLagrangeMultipliers(BSplinePatches{idI}.EFTPatches,lambda.EFTLagrangeMultipliers) + ...
        LambdaI;
    
    % Assemble the Lagrange multipliers field for the traction moments  of
    % the variational problem in patch I
    if isLagrangeMultipliersRotationsEnabled
        KLagrangeMultipliers(BSplinePatches{idI}.EFTPatches,mu.EFTLagrangeMultipliers) = ...
            KLagrangeMultipliers(BSplinePatches{idI}.EFTPatches,mu.EFTLagrangeMultipliers) + ...
            MuI;
    end
    
    % Assemble the Lagrange multipliers field for the traction forces  of
    % the variational problem in patch J
    KLagrangeMultipliers(BSplinePatches{idJ}.EFTPatches,lambda.EFTLagrangeMultipliers) = ...
        KLagrangeMultipliers(BSplinePatches{idJ}.EFTPatches,lambda.EFTLagrangeMultipliers) + ...
        LambdaJ;
    
    % Assemble the Lagrange multipliers field for the traction moments  of
    % the variational problem in patch J
    if isLagrangeMultipliersRotationsEnabled
        KLagrangeMultipliers(BSplinePatches{idJ}.EFTPatches,mu.EFTLagrangeMultipliers) = ...
            KLagrangeMultipliers(BSplinePatches{idJ}.EFTPatches,mu.EFTLagrangeMultipliers) + ...
            MuJ;
    end
    
    %% 1x. Assemble the Lagrange Multipliers matrices to master saddle point matrix corrssponding to the variational forms each Lagrange Multipliers field
    
    % Assemble the Lagrange Mulipliers matrix of patch I to the variational 
    % problem of the Lagrange Multipliers field of the traction forces
    KLagrangeMultipliers(lambda.EFTLagrangeMultipliers,BSplinePatches{idI}.EFTPatches) = ...
        KLagrangeMultipliers(lambda.EFTLagrangeMultipliers,BSplinePatches{idI}.EFTPatches) + ...
        LambdaI';
    
    % Assemble the Lagrange Mulipliers matrix of patch I to the variational 
    % problem of the Lagrange Multipliers field of the traction moments
    if isLagrangeMultipliersRotationsEnabled
        KLagrangeMultipliers(mu.EFTLagrangeMultipliers,BSplinePatches{idI}.EFTPatches) = ...
            KLagrangeMultipliers(mu.EFTLagrangeMultipliers,BSplinePatches{idI}.EFTPatches) + ...
            MuI';
    end
    
        % Assemble the Lagrange Mulipliers matrix of patch J to the variational 
    % problem of the Lagrange Multipliers field of the traction forces
    KLagrangeMultipliers(lambda.EFTLagrangeMultipliers,BSplinePatches{idJ}.EFTPatches) = ...
        KLagrangeMultipliers(lambda.EFTLagrangeMultipliers,BSplinePatches{idJ}.EFTPatches) + ...
        LambdaJ';
    
    % Assemble the Lagrange Mulipliers matrix of patch J to the variational
    % problem of the Lagrange Multipliers field of the traction moments
    if isLagrangeMultipliersRotationsEnabled
        KLagrangeMultipliers(mu.EFTLagrangeMultipliers,BSplinePatches{idJ}.EFTPatches) = ...
            KLagrangeMultipliers(mu.EFTLagrangeMultipliers,BSplinePatches{idJ}.EFTPatches) + ...
            MuJ';
    end
end

end
