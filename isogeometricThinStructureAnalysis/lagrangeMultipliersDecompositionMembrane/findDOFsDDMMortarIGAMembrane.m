%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _______________________________________________________               %
%   _______________________________________________________               %
%                                                                         %
%   Technische Universität München                                        %
%   Lehrstuhl für Statik, Prof. Dr.-Ing. Kai-Uwe Bletzinger               %
%   _______________________________________________________               %
%   _______________________________________________________               %
%                                                                         %
%                                                                         %
%   Authors                                                               %
%   _______________________________________________________________       %
%                                                                         %
%   Dipl.-Math. Andreas Apostolatos    (andreas.apostolatos@tum.de)       %
%   Dr.-Ing. Roland Wüchner            (wuechner@tum.de)                  %
%   Prof. Dr.-Ing. Kai-Uwe Bletzinger  (kub@tum.de)                       %
%   _______________________________________________________________       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [BSplinePatches,connections] = findDOFsDDMMortarIGAMembrane...
    (BSplinePatches,connections,propCoupling,fixedDOFs)
%% Function documentation
%
% Returns the enriched with the master and the slave DOFs for each patch
% pair connection structure for the enforcement of the multipatch coupling
% if membrane multipatches with the mortar method. Additionally each
% BSplinePatch array is enriched with the domain DOFs, that is, the DOFs
% which are not on any interface.
%
%          Input :
% BSplinePatches : Array of B-Spline patches each of which containing,
%                   .p,.q : The polynomial orders in each parametric
%                           direction
%                .Xi,.Eta : The knot vectors in each parametric direction
%                     .CP : The set of Control Point coordinates and
%                           weights
%    connections : Define the connection between the patches:
%                      .No : Number of connections
%               .xiEtaCoup : [patchID1 patchID2 xi12 eta12 xi21 eta21
%                             ...      ...    ...   ...  ...   ...]
%   propCoupling : Properties for the multipatch coupling using mortar
%                   .isSlaveSideCoarser : True if the coarser side of the
%                                         interface is chosen as slave,
%                                         otherwise the finer becomes the
%                                         slave
%      fixedDOFs : The global numbering of the DOFs which are constrained
%                  to a given value
%
%         Output :
%    connections : The enriched connections structure with
%                      .mortar : [idConnection idMaster idSlave]
%                  .masterDOFs : Structure array containing the master DOFs
%                                for the given connection
%                   .slaveDOFs : Structure array containing the slave DOFs
%                                for the given connection
% BSplinePatches : The enriched array of B-Spline patches with
%                  .domainDOFs : The global numbering in the multipatch 
%                                geometry of the DOFs which are not lying
%                                in any of the interface
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the connections
% ->
%    1i. Get the IDs of the patches involved
%
%   1ii. Create a basic Lagrange Multipliers field
%
%  1iii. Find the interface parametrization for the involved patches
%
%   1iv. Find the number of the interface DOFs from each patch
%
%    1v. Decide upon the master and the slave side
%
%   1vi. Find the patch indices in the multipatch geometry for both the master and the slave DOFs
%
%  1vii. Find the corresponding IDs of the interface DOFs
%
% 1viii. Degree elevate the Lagrange Multipliers field to match the slave side
%
%   1ix. Perform knot insertion to the Lagrange Multipliers field to match the slave side
%
%    1x. Add the Lagrange Multipliers field into the connections
% <-
%
% 2. Loop over all the B-Spline patches to get the purely domain DOFs
% ->
%    2i. Initialize the domain DOFs for the given patch
%
%   2ii. Find the connection into which the B-Spline patch belongs to
%
%  2iii. Loop over all the connections of the patch
%  ->
%        2iii.1. Get the coupling extension of the patch boundary for the given patch connection
%
%        2iii.2. Find the DOFs belonging to the given patch interface
%
%        2iii.3. Exclude the interface DOFs from the domain DOFs of the patch
%  <-
%
%   2iv. Assign the domain DOFs into the patch array
% <-
%
%% Function main body

%% 0. Read input

% Number of patches
noPatches = length(BSplinePatches);

% Initialize auxiliary arrays in case the mortar method is used and clear
% the Lagrange Multipliers field
connections.mortar = zeros(connections.No,2);
if isfield(connections,'lambda')
    connections = rmfield(connections,'lambda');
end
if isfield(connections,'mu')
    connections = rmfield(connections,'mu');
end
connections.masterDOFs = struct([]);
connections.slaveDOFs = struct([]);

% Clear the Lagrange Multipliers field if already defined
if isfield(connections,'lambda')
    connections = rmfield(connections,'lambda');
end

%% 1. Loop over all the connections
for iConnections = 1:connections.No
    %% 1i. Get the IDs of the patches involved
    idI = connections.xiEtaCoup(iConnections,1);
    idJ = connections.xiEtaCoup(iConnections,2);

    %% 1ii. Create a basic Lagrange Multipliers field
    pLambda = 0;
    XiLambda = [0 1];
    CPLambda(:,4) = 1;
    isNURBSLambda = 0;
    nxiLambda = length(CPLambda(:,1,1));
    for i = 1:nxiLambda
        if CPLambda(i,4)~=1
            isNURBSLambda = 1;
            break;
        end
    end

    %% 1iii. Find the interface parametrization for the involved patches
    xiCoupI = connections.xiEtaCoup(iConnections,3:4);
    etaCoupI = connections.xiEtaCoup(iConnections,5:6);
    if xiCoupI(1,1) ~= xiCoupI(1,2) && etaCoupI(1,1) == etaCoupI(1,2)
        isOnXiI = true;
    elseif xiCoupI(1,1) == xiCoupI(1,2) && etaCoupI(1,1) ~= etaCoupI(1,2)
        isOnXiI = false;
    else
        error('Either the xi or the eta direction has to be fixed along the interface for patch %d',idI);
    end
    xiCoupJ = connections.xiEtaCoup(iConnections,7:8);
    etaCoupJ = connections.xiEtaCoup(iConnections,9:10);
    if xiCoupJ(1,1) ~= xiCoupJ(1,2) && etaCoupJ(1,1) == etaCoupJ(1,2)
        isOnXiJ = true;
    elseif xiCoupJ(1,1) == xiCoupJ(1,2) && etaCoupJ(1,1) ~= etaCoupJ(1,2)
        isOnXiJ = false;
    else
        error('Either the xi or the eta direction has to be fixed along the interface for patch %d',idJ);
    end
    
    %% 1iv. Find the number of the interface DOFs from each patch
%     if isOnXiI
%         noCPsI = length(BSplinePatches{idI}.CP(:,1,1));
%     else
%         noCPsI = length(BSplinePatches{idI}.CP(1,:,1));
%     end
%     noDOFsI = 3*noCPsI;
%     if isOnXiJ
%         noCPsJ = length(BSplinePatches{idJ}.CP(:,1,1));
%     else
%         noCPsJ = length(BSplinePatches{idJ}.CP(1,:,1));
%     end
%     noDOFsJ = 3*noCPsJ;
    
    %% 1v. Decide upon the master and the slave side
%     isFlipped = false;
%     cond = noDOFsI < noDOFsJ;
%     if ~propCoupling.isSlaveSideCoarser
%         cond = noDOFsI >= noDOFsJ;
%     end         
%     if cond
%         isFlipped = true;
%     end
    
    isFlipped = ~connections.masterSlave(iConnections,1);

    connections.mortar(iConnections,:) = [idI idJ];
    if isFlipped
        connections.mortar(iConnections,:) = fliplr(connections.mortar(iConnections,:));
    end
    
    %% 1vi. Find the patch indices in the multipatch geometry for both the master and the slave DOFs
    indexI = BSplinePatches{idI}.EFTPatches(1,1) - 1;
    indexJ = BSplinePatches{idJ}.EFTPatches(1,1) - 1;
    indexMaster = indexI;
    indexSlave = indexJ;
    if isFlipped
        indexMaster = indexJ;
        indexSlave = indexI;
    end
    
    %% 1vii. Find the corresponding IDs of the interface DOFs
    cbI = [];
    for dir = 1:3
        cbI = findDofs3D(cbI,xiCoupI,etaCoupI,dir,BSplinePatches{idI}.CP);
    end
    cbJ = [];
    for dir = 1:3
        cbJ = findDofs3D(cbJ,xiCoupJ,etaCoupJ,dir,BSplinePatches{idJ}.CP);
    end
    cbMaster = cbI;
    cbSlave = cbJ;
    if isFlipped
        cbMaster = cbJ;
        cbSlave = cbI;
    end
    connections.masterDOFs{iConnections} = indexMaster + cbMaster;
    connections.slaveDOFs{iConnections} = indexSlave + cbSlave;
    
    %% 1viii. Degree elevate the Lagrange Multipliers field to match the slave side
    isOnXi = isOnXiJ;
    if isFlipped
        isOnXi = isOnXiI;
    end
    pLM = BSplinePatches{connections.mortar(iConnections,2)}.q;
    if isOnXi
        pLM = BSplinePatches{connections.mortar(iConnections,2)}.p;
    end
    if pLM > 0
        clear pLambda XiLambda CPLambda;
        pLambda = 1;
        XiLambda = [0 0 1 1];
        CPLambda(:,4) = [1 1];
        isNURBSLambda = 0;
        nxiLambda = length(CPLambda(:,1,1));
        for i = 1:nxiLambda
            if CPLambda(i,4) ~= 1
                isNURBSLambda = 1;
                break;
            end
        end
        
        % Perform accordingly a p-refinement
        tpLambda = pLM - pLambda;
        [XiLambda,CPLambda,pLambda] = degreeElevateBSplineCurve...
            (pLambda,XiLambda,CPLambda,tpLambda,'');
    end
    
    %% 1ix. Perform knot insertion to the Lagrange Multipliers field to match the slave side
    interiorKnots = BSplinePatches{connections.mortar(iConnections,2)}.Eta;
    if isOnXi
        interiorKnots = BSplinePatches{connections.mortar(iConnections,2)}.Xi;
    end
    noKnots = length(interiorKnots);
    interiorKnots = interiorKnots(pLM + 2:noKnots - pLM - 1);
    [XiLambdaLM,CPLambdaLM] = knotRefineBSplineCurve...
        (pLambda,XiLambda,CPLambda,interiorKnots,'');
    
    %% 1x. Add the Lagrange Multipliers field into the connections
    connections.lambda{iConnections} = fillUpLagrangeMultipliers...
        (pLambda,XiLambdaLM,CPLambdaLM,isNURBSLambda);
end

%% 2. Loop over all the B-Spline patches to get the purely domain DOFs
% for iPatches = 1:noPatches
%     %% 2i. Initialize the domain DOFs for the given patch
%     domainDOFs = BSplinePatches{iPatches}.EFTPatches;
%     
%     %% 2ii. Find the connection into which the B-Spline patch belongs to
%     [indexI,indexJ] = find(iPatches == connections.mortar);
%     if length(indexI) ~= length(indexJ)
%         error('Indices of patch connections are not matching');
%     end
%     noPatchConnections = length(indexI);
%     
%     %% 2iii. Loop over all the connections of the patch
%     for iPatchConn = 1:noPatchConnections
%         %% 2iii.1. Get the coupling extension of the patch boundary for the given patch connection
%         idConn = indexI(iPatchConn,1);
%         idPos = indexJ(iPatchConn,1);
%         if idPos ~= 1 && idPos ~= 2
%             error('The position of the patch id must be either 1 or 2');
%         end
%         isMaster = true;
%         if idPos == 2
%             isMaster = false;
%         end
%         
%         %% 2iii.2. Find the DOFs belonging to the given patch interface
%         interfaceDOFs = connections.masterDOFs{idConn};
%         if ~isMaster
%             interfaceDOFs = connections.slaveDOFs{idConn};
%         end
%         
%         %% 2iii.3. Exclude the interface DOFs from the domain DOFs of the patch
%         domainDOFs = domainDOFs(~ismember(domainDOFs,interfaceDOFs));
%     end
%     
%     %% 2iv. Assign the domain DOFs into the patch array
%     BSplinePatches{iPatches}.domainDOFs = domainDOFs;
% end

%% 3. Clean the master and the slave DOF arrays from the DOFs where homogeneous Dirichlet conditions are applied
for iConnOuter = 1:connections.No
    % Apply the homogeneous Dirichlet boundary conditions at the master
    % DOFs
    connections.masterDOFs{iConnOuter} = connections.masterDOFs{iConnOuter}(~ismember(connections.masterDOFs{iConnOuter},fixedDOFs));

    % Find the slave DOFs which are located on the Dirichlet boundary
    slaveHomDOFs = ismember(connections.slaveDOFs{iConnOuter},fixedDOFs);

    % Apply the homogeneous Dirichlet boundary conditions at the slave
    % DOFs
    connections.slaveDOFs{iConnOuter} = connections.slaveDOFs{iConnOuter}(~slaveHomDOFs);

    %% Loop over all previous connections and check if the slave DOFs have been aready defined in another connection
%     for iConnInner = 1:iConnOuter - 1
%         % Check if there are slave DOFs that are defined as master DOFs
%         % in previous connections and exclude them
%         slaveMasterDOFs = ismember(connections.slaveDOFs{iConnOuter},connections.masterDOFs{iConnInner});
%         connections.slaveDOFs{iConnOuter} = connections.slaveDOFs{iConnOuter}(~slaveMasterDOFs);
% 
%         % Check if there are master DOFs that are defined as slave DOFs
%         % in previous connections and exclude them
%         masterSlaveDOFs = ismember(connections.masterDOFs{iConnOuter},connections.slaveDOFs{iConnInner});
%         connections.slaveDOFs{iConnInner} = connections.slaveDOFs{iConnInner}(~masterSlaveDOFs);
% 
% 
% % %             slaveCommonDOFs = ismember(connections.slaveDOFs{iConnOuter},connections.slaveDOFs{iConnInner});
% %             slaveCommonDOFs = ismember(connections.slaveDOFs{iConnOuter},connections.masterDOFs{iConnInner});
% %             
% %             connections.slaveDOFs{iConnOuter} = connections.slaveDOFs{iConnOuter}(~slaveCommonDOFs);
%     end
end

%% 4. Clean the slave DOF arrays from DOFs which are defined as slave DOFs in one connection but as master DOFs in another
for iConnOuter = 1:connections.No
    for iConnInner = 1:connections.No 
%         % Check if there are slave DOFs that are defined as master DOFs
%         % in previous connections and exclude them
%         slaveMasterDOFs = ismember(connections.slaveDOFs{iConnOuter},connections.masterDOFs{iConnInner});
%         connections.slaveDOFs{iConnOuter} = connections.slaveDOFs{iConnOuter}(~slaveMasterDOFs);

        masterSlaveDOFs = ismember(connections.masterDOFs{iConnOuter},connections.slaveDOFs{iConnInner});
        connections.masterDOFs{iConnOuter} = connections.masterDOFs{iConnOuter}(~masterSlaveDOFs);
    end
end

end