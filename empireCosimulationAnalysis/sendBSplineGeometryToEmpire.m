function sendBSplineGeometryToEmpire(BSplinePatches,connections,...
    propEmpireCoSimulation,tab,outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Sends a multipatch B-Spline surface to Empire for co-simulation.
%
%                  Input :
%         BSplinePatches : Structure containing all the information 
%                          regarding the onnections between the multipatches
%            connections : Define the connection between the patches:
%                            .No : Number of connections
%                     .xiEtaCoup : [patchID1 patchID2 xi12 eta12 xi21 eta21
%                                   ...      ...    ...   ...  ...   ...]
% propEmpireCoSimulation : Properties for the co-simulation with EMPIRE
%                            .isCoSimulation : Flag on whether 
%                                              co-simulation with EMPIRE is 
%                                              assumed
%                          .isInterfaceLayer : Flag on whether the matlab 
%                                              client is used as an 
%                                              interface layer
%             .isGaussPointProvidedDirichlet : Flag on whether the Gauss
%                                              Point data along the
%                                              Dirichlet boundaries are
%                                              provided
%             .isGaussPointProvidedInterface : Flag on whether the Gauss
%                                              Point data along the
%                                              interface boundaries are
%                                              provided
%                          .propIntDirichlet : Integration along the
%                                              Dirichlet boundaries
%                                              properties,
%                                                .type : 'default', 'user'
%                                               .noGPs : Number of Gauss
%                                                        Points
%                          .propIntInterface : Integration along the
%                                              interface boundaries
%                                              properties,
%                                                .type : 'default', 'user'
%                                               .noGPs : Number of Gauss
%                                                        Points
%                              .strMatlabXml : Name of the xml file for 
%                                              connection to Empire
%                    tab : Tabulation when printing message in the command 
%                          window
%                 outMsg : Enables message outputting in the command window 
%                          if chosen as 'outputEnabled'
%
%         Output :
%                  Sends B-Splint geometry to Empire
%
% Function layout :
%
% 0. Read input
%
% 1. Get the total number of Control Points in the multipatch geometry
%
% 2. Send the mesh type to Empire
%
% 3. Loop over all patches
% ->
%    3i. Get patch parameters
%
%   3ii. Initialize arrays containing the numbering of the Control Points and the unique numbering of Control Points for each patch
%
%  3iii. Loop over all Control Points
%  ->
%        3iii.1. Update outer counter
%
%        3iii.2. Assign the Control Point coordinates and weights to the Control Point net array
%
%        3iii.3. Update the inner counter
%
%        3iii.4. Create the array of the unique IDs starting from zero
%  <-
%
%   3iv. Send the patch to Empire
%
%    3v. Send the corresponding trimming information to Empire
%
%   3vi. Loop over all the trimming loops of the patch
%   ->
%        3vi.1. Send the trimming loop info
%
%        3vi.2. Loop over all trimming curves of the trimming loop and send the trimming gometrical information
%   <-
% <-
%
% 4. Send the number of Dirichlet boundary conditions to Empire
%
% 5. Loop over all patches
% ->
%    5i. Get patch counter
%
%   5ii. Loop over all Dirichlet boundaries of the patch
%   ->
%        5ii.1. Get the extensions of the Dirichlet boundary
%
%        5ii.2. Send the description of the Dirichlet boundary along which weak imposition of constraints is assumed
%
%        5ii.3. Compute the coupling data along the Dirichlet boundary where conditions are weakly applied
%
%        5ii.4. Send the coupling data along the Dirichlet boundary where conditions are weakly applied
%   <-
% <-
%
% 6. Send the number of patch connections to Empire
%
% 7. Loop over all connections
% ->
%    7i. Get the counters of the master and slave patches
%
%   7ii. Get the extensions of the interface boundary for the master patch
%
%  7iii. Get the extensions of the interface boundary for the slave patch
%
%   7iv. Send the description of the interface boundary along which weak continuity conditions are applied
%
%    7v. Check the orientation of the interface parametrizations from both patches
%
%   7vi. Compute the coupling data along the interface boundary where conditions are weakly applied
%
%  7vii. Send the coupling data along the interface boundary where conditions are weakly applied
% <-
%
% 8. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf(['\n' tab,'Sending geometrical information to Empire\n\n']);
end

%% 0. Read input

% Flag on whether the Gauss Point data for the Dirichlet boundary
% conditions are provided
assert(isfield(propEmpireCoSimulation,'isGaussPointProvidedDirichlet'));
if propEmpireCoSimulation.isGaussPointProvidedDirichlet
    isDirichletGPProvided = 1;
    assert(isfield(propEmpireCoSimulation,'propIntDirichlet'));
    assert(isfield(propEmpireCoSimulation.propIntDirichlet,'type'));
    assert(isfield(propEmpireCoSimulation.propIntDirichlet,'noGPs'));
    propIntDirichlet = propEmpireCoSimulation.propIntDirichlet;
else
    isDirichletGPProvided = 0;
end

% Flag on whether the Gauss Point data for the interface continuity
% conditions are provided
assert(isfield(propEmpireCoSimulation,'isGaussPointProvidedInterface'));
if propEmpireCoSimulation.isGaussPointProvidedInterface
    isCouplingGPProvided = 1;
    assert(isfield(propEmpireCoSimulation,'propIntInterface'));
    assert(isfield(propEmpireCoSimulation.propIntInterface,'type'));
    assert(isfield(propEmpireCoSimulation.propIntInterface,'noGPs'));
    propIntInterface = propEmpireCoSimulation.propIntInterface;
else
    isCouplingGPProvided = 0;
end

% Get the properties of the Dirichlet and interface integration arrays

% Name of the mesh to be sent
meshName = 'IGAMesh';

% Number of Patches
noPatches = length(BSplinePatches);

% Flag on the trimming
isTrimmed = 1;

% Flag on whether the trimming loop is an inner or outer loop
isInner = 0;

% Counter for the boundary loop of each patch
patchBLCtr = 0;

% Number of boundary loops per patch
noLoopsPatch = 1;

% Number of trimming curves per patch
noTrimmingCurvesPatch = 4;

% Geometrical information of the trimming curve for each patch
dirTC = 1;
pTC = 1;
noKnotsTC = 4;
knotVctTC = [0 0 1 1];
noCPsTC = 2;
%
%      eta
%       ^
%       |
%       |     (2)
%       -------<-------
%       |             |
%       |             |
%  (3)  v             ^ (1)
%       |             |
%       |             |
%       ------->------------> xi
%             (0)
%
cpNetTC = [0 0 0 1 1 0 0 1
           1 0 0 1 1 1 0 1
           1 1 0 1 0 1 0 1
           0 1 0 1 0 0 0 1];

% Initialize counter
counterOuter = 0;

%% 1. Get the total number of Control Points in the multipatch geometry
noCPs = 0;
for iPatches = 1:noPatches
  noCPs = noCPs + BSplinePatches{iPatches}.noCPs;
end

%% 2. Send the mesh type to Empire
if strcmp(outMsg,'outputEnabled')
    fprintf(['\n' tab,'Sending %d patches to Empire\n\n'],noPatches);
end
EMPIRE_API_sendIGAMesh(meshName,noPatches,noCPs);

%% 3. Loop over all patches
for iPatches = 1:noPatches 
    %% 3i. Get patch parameters
    p = BSplinePatches{iPatches}.p;
    q = BSplinePatches{iPatches}.q;
    noKnotsXi = length(BSplinePatches{iPatches}.Xi);
    noKnotsEta = length(BSplinePatches{iPatches}.Eta);
    Xi = BSplinePatches{iPatches}.Xi;
    Eta = BSplinePatches{iPatches}.Eta;
    noXiCP = length(BSplinePatches{iPatches}.CP(:,1,1));
    noEtaCP = length(BSplinePatches{iPatches}.CP(1,:,1));

    %% 3ii. Initialize arrays containing the numbering of the Control Points and the unique numbering of Control Points for each patch
    netCP = zeros(1,4*BSplinePatches{iPatches}.noCPs);
    netUniqueIDs = zeros(1,BSplinePatches{iPatches}.noCPs);
    counterInner = 0;
    
    %% 3iii. Loop over all Control Points
    for j = 1:noEtaCP
        for k = 1:noXiCP
            %% 3iii.1. Update outer counter
            counterOuter = counterOuter + 1;

            %% 3iii.2. Assign the Control Point coordinates and weights to the Control Point net array
            netCP(1,4*counterInner + 1) = BSplinePatches{iPatches}.CP(k,j,1);
            netCP(1,4*counterInner + 2) = BSplinePatches{iPatches}.CP(k,j,2);
            netCP(1,4*counterInner + 3) = BSplinePatches{iPatches}.CP(k,j,3);
            netCP(1,4*counterInner + 4) = BSplinePatches{iPatches}.CP(k,j,4);

            %% 3iii.3. Update the inner counter
            counterInner = counterInner + 1;
            
            %% 3iii.4. Create the array of the unique IDs starting from zero
            netUniqueIDs(1,counterInner) = counterOuter - 1;
        end        
    end

    %% 3iv. Send the patch to Empire
    if strcmp(outMsg,'outputEnabled')
        fprintf([tab '\t' 'Sending patch %d/%d to Empire\n'],iPatches,noPatches);
    end
    EMPIRE_API_sendIGAPatch(p,noKnotsXi,Xi,q,noKnotsEta,Eta,noXiCP,noEtaCP,netCP,netUniqueIDs);

    %% 3v. Send the corresponding trimming information to Empire
    if strcmp(outMsg,'outputEnabled')
        fprintf([tab '\t' '\t' 'Sending trimming information to Empire\n']);
    end
    EMPIRE_API_sendIGATrimmingInfo(isTrimmed,noLoopsPatch);
    
    %% 3vi. Loop over all the trimming loops of the patch
    for iTrimmingLoop = 1:noLoopsPatch
        %% 3vi.1. Send the trimming loop info
        EMPIRE_API_sendIGATrimmingLoopInfo(isInner,noTrimmingCurvesPatch);
        
        %% 3vi.2. Loop over all trimming curves of the trimming loop and send the trimming gometrical information
        for iTrimmingCurves = 1:noTrimmingCurvesPatch
            if strcmp(outMsg,'outputEnabled')
                fprintf([tab '\t' '\t' 'Sending trimming curve %d/%d to Empire\n'],iTrimmingCurves,noTrimmingCurvesPatch);
            end
            EMPIRE_API_sendIGATrimmingCurve(dirTC,pTC,noKnotsTC,knotVctTC,noCPsTC,cpNetTC(iTrimmingCurves,:));
        end
    end
end

%% 4. Send the number of Dirichlet boundary conditions to Empire
noWeakDBC = 0;
for iPatches = 1:noPatches
    noWeakDBC = noWeakDBC + BSplinePatches{iPatches}.weakDBC.noCnd;
end
if strcmp(outMsg,'outputEnabled')
    fprintf(['\n' tab 'Sending %d Dirichlet conditions to Empire\n'],noWeakDBC);
end
EMPIRE_API_sendIGANumDirichletConditions(noWeakDBC);

%% 5. Loop over all patches
for iPatches = 1:noPatches
    %% 5i. Get patch counter
    patchCtr = iPatches - 1;
    
    %% 5ii. Loop over all Dirichlet boundaries of the patch
    for iCnd = 1:BSplinePatches{iPatches}.weakDBC.noCnd
        %% 5ii.1. Get the extensions of the Dirichlet boundary
        xiExtension = BSplinePatches{iPatches}.weakDBC.xiExtension{iCnd};
        etaExtension = BSplinePatches{iPatches}.weakDBC.etaExtension{iCnd};
        if etaExtension(1,1) == etaExtension(1,2)
            if etaExtension(1,1) == 0
                patchBLTrCurveCtr = 0;
            elseif etaExtension(1,1) == 1
                patchBLTrCurveCtr = 2;
            else
                error('The extension along eta for the Dirichlet boundary must be either 0 or 1');
            end
        elseif xiExtension(1,1) == xiExtension(1,2)
            if xiExtension(1,1) == 1
                patchBLTrCurveCtr = 1;
            elseif xiExtension(1,1) == 0
                patchBLTrCurveCtr = 3;
            else
                error('The extension along xi for the Dirichlet boundary must be either 0 or 1');
            end
        else
           error('One of the extensions, xi or eta, has to be fixed for the description of a Dirichlet boundary'); 
        end
        
        %% 5ii.2. Send the description of the Dirichlet boundary along which weak imposition of constraints is assumed
        if strcmp(outMsg,'outputEnabled')
            fprintf([tab '\t' 'Sending trimming curve info along Dirichlet boundary [%d,%d]x[%d,%d] of patch %d to Empire\n'],xiExtension(1,1),xiExtension(1,2),etaExtension(1,1),etaExtension(1,2),iPatches);
        end
        EMPIRE_API_sendIGADirichletConditionInfo(patchCtr,patchBLCtr,patchBLTrCurveCtr,isDirichletGPProvided);
        
        %% 5ii.3. Compute the coupling data along the Dirichlet boundary where conditions are weakly applied
        if isDirichletGPProvided
            [trCurveNumGP, trCurveGPs, trCurveGPWeights, trCurveGPTangents, trCurveGPJacobianProducts] = ...
                computeGaussPointDataAlongDirichletBoundary...
                (BSplinePatches{iPatches},xiExtension,etaExtension,propIntDirichlet,tab,outMsg);
        end
        
        %% 5ii.4. Send the coupling data along the Dirichlet boundary where conditions are weakly applied
        if isDirichletGPProvided
            EMPIRE_API_sendIGADirichletConditionData(trCurveNumGP, trCurveGPs, trCurveGPWeights, trCurveGPTangents, trCurveGPJacobianProducts);
        end
    end
end

%% 6. Send the number of patch connections to Empire
isConnections = false;
if strcmp(outMsg,'outputEnabled')
    if ~ischar(connections)
        if ~isempty(connections)
            if isfield(connections,'No')
                isConnections = true;
            end
        end
    end
end
if isConnections
    noConnections = connections.No;
else
    noConnections = 0;
end
fprintf(['\n' tab 'Sending %d patch connections to Empire\n'],noConnections);
EMPIRE_API_sendIGANumPatchConnections(noConnections);

%% 7. Loop over all connections
for iConn = 1:noConnections
    %% 7i. Get the counters of the master and slave patches
    IDMaster = connections.xiEtaCoup(iConn,1);
    masterPatchCtr = IDMaster - 1;
    IDSlave = connections.xiEtaCoup(iConn,2);
    slavePatchCtr = IDSlave - 1;
    
    %% 7ii. Get the extensions of the interface boundary for the master patch
    xiExtensionMaster = connections.xiEtaCoup(iConn,3:4);
    etaExtensionMaster = connections.xiEtaCoup(iConn,5:6);
    if etaExtensionMaster(1,1) == etaExtensionMaster(1,2)
        if etaExtensionMaster(1,1) == 0
            masterPatchBLTrCurveCtr = 0;
        elseif etaExtensionMaster(1,1) == 1
            masterPatchBLTrCurveCtr = 2;
        else
            error('The extension along eta for the Dirichlet boundary must be either 0 or 1');
        end
    elseif xiExtensionMaster(1,1) == xiExtensionMaster(1,2)
        if xiExtensionMaster(1,1) == 1
            masterPatchBLTrCurveCtr = 1;
        elseif xiExtensionMaster(1,1) == 0
            masterPatchBLTrCurveCtr = 3;
        else
            error('The extension along xi for the Dirichlet boundary must be either 0 or 1');
        end
    else
       error('One of the extensions, xi or eta, has to be fixed for the description of a Dirichlet boundary'); 
    end
    
    %% 7iii. Get the extensions of the interface boundary for the slave patch
    xiExtensionSlave = connections.xiEtaCoup(iConn,7:8);
    etaExtensionSlave = connections.xiEtaCoup(iConn,9:10);
    if etaExtensionSlave(1,1) == etaExtensionSlave(1,2)
        if etaExtensionSlave(1,1) == 0
            slavePatchBLTrCurveCtr = 0;
        elseif etaExtensionSlave(1,1) == 1
            slavePatchBLTrCurveCtr = 2;
        else
            error('The extension along eta for the Dirichlet boundary must be either 0 or 1');
        end
    elseif xiExtensionSlave(1,1) == xiExtensionSlave(1,2)
        if xiExtensionSlave(1,1) == 1
            slavePatchBLTrCurveCtr = 1;
        elseif xiExtensionSlave(1,1) == 0
            slavePatchBLTrCurveCtr = 3;
        else
            error('The extension along xi for the Dirichlet boundary must be either 0 or 1');
        end
    else
       error('One of the extensions, xi or eta, has to be fixed for the description of a Dirichlet boundary'); 
    end
    
    %% 7iv. Send the description of the interface boundary along which weak continuity conditions are applied
    if strcmp(outMsg,'outputEnabled')
        fprintf([tab '\t' 'Sending trimming curve info along interface boundary [%d,%d]x[%d,%d] of patch %d with [%d,%d]x[%d,%d] of patch %d to Empire\n'],xiExtensionMaster(1,1),xiExtensionMaster(1,2),etaExtensionMaster(1,1),etaExtensionMaster(1,2),masterPatchCtr + 1,...
            xiExtensionSlave(1,1),xiExtensionSlave(1,2),etaExtensionSlave(1,1),etaExtensionSlave(1,2),slavePatchCtr + 1);
    end
    EMPIRE_API_sendIGAPatchConnectionInfo(masterPatchCtr,patchBLCtr,masterPatchBLTrCurveCtr,...
                                          slavePatchCtr,patchBLCtr,slavePatchBLTrCurveCtr,...
                                          isCouplingGPProvided);
                                      
    %% 7v. Check the orientation of the interface parametrizations from both patches
    if isCouplingGPProvided
        haveSameOrientation = findSubdomainInterfaceOrientation...
            (BSplinePatches{IDMaster}.p,BSplinePatches{IDMaster}.Xi,BSplinePatches{IDMaster}.q,BSplinePatches{IDMaster}.Eta,BSplinePatches{IDMaster}.CP,BSplinePatches{IDMaster}.isNURBS,xiExtensionMaster,etaExtensionMaster,...
            BSplinePatches{IDSlave}.p,BSplinePatches{IDSlave}.Xi,BSplinePatches{IDSlave}.q,BSplinePatches{IDSlave}.Eta,BSplinePatches{IDSlave}.CP,BSplinePatches{IDSlave}.isNURBS,xiExtensionSlave,etaExtensionSlave);
    end
                                      
    %% 7vi. Compute the coupling data along the interface boundary where conditions are weakly applied
    if isCouplingGPProvided
        [couplingCurveNumGP, trCurveMasterGPs, trCurveSlaveGPs, trCurveGPWeights, trCurveMasterGPTangents, trCurveSlaveGPTangents, trCurveGPJacobianProducts] = ...
            computeGaussPointDataAlongInterfaceBoundary...
            (BSplinePatches{connections.xiEtaCoup(iConn,1)},BSplinePatches{connections.xiEtaCoup(iConn,2)},xiExtensionMaster,etaExtensionMaster,xiExtensionSlave,etaExtensionSlave,haveSameOrientation,propIntInterface,tab,outMsg);
    end
    
    %% 7vii. Send the coupling data along the interface boundary where conditions are weakly applied
    if isCouplingGPProvided
        EMPIRE_API_sendIGAPatchConnectionData(couplingCurveNumGP, trCurveMasterGPs, trCurveSlaveGPs, trCurveGPWeights, trCurveMasterGPTangents, trCurveSlaveGPTangents, trCurveGPJacobianProducts);
    end
end

%% 8. Appendix
if strcmp(outMsg,'outputEnabled')
    fprintf('\n');
end

end
