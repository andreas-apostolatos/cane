function TMortar = computeConstantMtxForDDMMortarIGAThinStructure ...
    (BSplinePatches, connections, numDOFs, propCoupling)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the mortar transformation matrices for each patch pair in a cell 
% array corresponding to the mortar method for the multipatch coupling of 
% thin-walled isogeometric patches.
%
%       BSplinePatches : Structure containing all the information regarding 
%                        the connections between the multipatches
%          connections : Define the connection between the patches:
%                            .No : Number of connections
%                     .xiEtaCoup : [patchID1 patchID2 xi12 eta12 xi21 eta21
%                                 ...      ...    ...   ...  ...   ...]
%               noDOFs : The complete number of DOFs
%         propCoupling : Properties of the multipatch coupling
%                             .intC : On the integration of the coupling
%                                     interface
%
%               Output :
%              TMortar : Cell array of the mortar transformation matrices
%                        at each patch pair interface
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the connections in the multipatch geometry
% ->
%    1i. Get the patch IDs
%
%   1ii. Get the coupling boundaries for both patches
%
%  1iii. Determine the interface orientation
%
%   1iv. Compute the Lagrange Multipliers matrices
%
%    1v. Compute the mortar transformation matrix for the given patch connection
% <-
%
%% Function main body

%% 0. Read input

% Initialize output array
TMortar = struct([]);

%% 1. Loop over all the connections in the multipatch geometry
for iConnections = 1:connections.No
    %% 1i. Get the patch IDs
    
    % Patch I :
    % _________
    
    idI = connections.xiEtaCoup(iConnections,1);
    
    % Patch J :
    % _________
    
    idJ = connections.xiEtaCoup(iConnections,2);
    
    %% 1ii. Get the coupling boundaries for both patches
    
    % Patch I :
    % _________
    
    BSplinePatches{idI}.xicoup = connections.xiEtaCoup(iConnections,3:4);
    BSplinePatches{idI}.etacoup = connections.xiEtaCoup(iConnections,5:6);
    
	% Patch J :
    % _________
    
    BSplinePatches{idJ}.xicoup = connections.xiEtaCoup(iConnections,7:8);
    BSplinePatches{idJ}.etacoup = connections.xiEtaCoup(iConnections,9:10);
    
    %% 1iii. Determine the interface orientation
    isSameOrientation = findSubdomainInterfaceOrientation...
        (BSplinePatches{idI}.p, BSplinePatches{idI}.Xi, BSplinePatches{idI}.q, ...
        BSplinePatches{idI}.Eta, BSplinePatches{idI}.CP, BSplinePatches{idI}.isNURBS, ...
        BSplinePatches{idI}.xicoup, BSplinePatches{idI}.etacoup, BSplinePatches{idJ}.p, ...
        BSplinePatches{idJ}.Xi, BSplinePatches{idJ}.q, BSplinePatches{idJ}.Eta, ...
        BSplinePatches{idJ}.CP, BSplinePatches{idJ}.isNURBS, BSplinePatches{idJ}.xicoup, ...
        BSplinePatches{idJ}.etacoup);
    
    %% 1iv. Compute the Lagrange Multipliers matrices
    mortarConnection = connections.mortar(iConnections,:);
    if idI == mortarConnection(1,1) && idJ == mortarConnection(1,2)
        isFlipped = false;
    elseif idI == mortarConnection(1,2) && idJ == mortarConnection(1,1)
        isFlipped = true;
    else
        error('The IDs of the master and slave patch do not match with the IDs stored in the mortar array');
    end
    indexI = BSplinePatches{idI}.EFTPatches(1,1) - 1;
    indexJ = BSplinePatches{idJ}.EFTPatches(1,1) - 1;
    indexMaster = indexI;
    indexSlave = indexJ;
    patchMaster = BSplinePatches{idI};
    patchSlave = BSplinePatches{idJ};
    if isFlipped
        indexMaster = indexJ;
        indexSlave = indexI;
        patchMaster = BSplinePatches{idJ};
        patchSlave = BSplinePatches{idI};
    end
    masterDOFs = connections.masterDOFs{iConnections} - indexMaster;
    slaveDOFs = connections.slaveDOFs{iConnections} - indexSlave;
    [LambdaMaster, LambdaSlave] = computeDDMLagrangeMultipliersMtces4MortarIGAThinStructure ...
        (patchMaster, patchSlave, isSameOrientation, propCoupling);
    LambdaMaster = LambdaMaster(masterDOFs,slaveDOFs)';
    LambdaSlave = LambdaSlave(slaveDOFs,slaveDOFs)';
    
    %% 1v. Compute the mortar transformation matrix for the given patch connection
    TMortar{iConnections} = - LambdaSlave\LambdaMaster;
end

end
