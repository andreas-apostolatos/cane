function F = computeLoadVctFEMThermalConductionAnalysis...
    (strMsh,analysis,propNBC,t,gaussInt,outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%                  Marko Leskovar
%
%% Function documentation
%
% Returns the global flux (load) vector corresponding to classical Finite 
% Element analysis for the case of the heat conduction problem.
%
%           Input :
%          strMsh : Nodes and elements in the mesh
%        analysis : .type : The analysis type
%         propNBC :    .nodes : The nodes where Neumann boundary 
%                                conditions are applied
%                    .loadType : The type of the load for each Neumann node
%                   .fctHandle : The function handle for each Neumann node
%                                for the computation of the load vector 
%                                (these functions are unde the folder load)
%               t : The time instance
%        gaussInt : On the numerical quadrature scheme
%                       int.domain :
%                         .type : 'default' or 'manual'
%                 .boundaryNoGP : The number of Gauss points for the
%                                 boundary integration
%          outMsg : On outputting information
%
%          Output :
%               F : The global load vector
%
% Function layout :
%
% 0. Read input
%
% 1. Get the Gauss integration rules
%
% 2. Loop over all elements on the Neumann boundary
% ->
%    2i. Get the nodes which are on the Neumann boundary
%
%   2ii. Get the nodes of the element on the Neumann boundary
%
%  2iii. Get the function handle for this type of loading
%
%   2iv. Assemble the Element Freedome Table (EFT)
%
%    2v. Loop over all Gauss Points
%    ->
%        2v.1. Transform the Gauss Point from the parameter to the physical space
%
%        2v.2. Compute the basis functions at the Gauss Point
%
%        2v.3. Sort the basis functions into the basis functions array
%
%        2v.4. Compute the applied load onto the Gauss Point
%
%        2v.5. Compute the determinant of the transformation from the physical space to the parameter space
%
%        2v.6. Compute the element load vector on the Gauss point
%
%        2v.7. Assemble the contribution to the global load vector via the EFT
%    <-
% <-
% 3. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('______________________________________________________________\n');
    fprintf('##############################################################\n');
    fprintf('Computation the global load vector for a FEM discretized plate\n');
    fprintf('in heat transfer problem has been initiated\n');
    fprintf('______________________________________________________________\n\n');
    tic;
end

%% 0. Read input

% Number of nodes
numNodes = length(strMsh.nodes(:, 1));

% Number of DOFs
if strcmp(analysis.type,'THERMAL_CONDUCTION_2D')
    nmuDOFs = numNodes;
else
    error('The analysis type must be THERMAL_CONDUCTION_2D');
end

% Initialize output array
F = zeros(nmuDOFs, 1);

%% 1. Get the Gauss integration rules

% For the boundary integration
if strcmp(gaussInt.type, 'default')
    numGP = 1;
elseif strcmp(gaussInt.type, 'user')
    numGP = gaussInt.boundaryNoGP;
end
[GP, GW] = getGaussPointsAndWeightsOverUnitDomain(numGP);

%% 2. Loop over all elements on the Neumann boundary
for iElmnt = 1:length(propNBC.lines(:, 1))
    %% 2i. Get the nodes which are on the Neumann boundary
    nodeIDs = propNBC.lines(iElmnt, 1:2);
    x1 = strMsh.nodes(nodeIDs(1), :);
    x2 = strMsh.nodes(nodeIDs(2), :);
    
    %% 2ii. Get the nodes of the element on the Neumann boundary
    elementID = propNBC.lines(iElmnt, 3);
    element = strMsh.elements(elementID, :);
    index = isnan(element);
    element(index) = [];
    numNodesEl = length(element);
    numDOFsEl = numNodesEl;
    nodes = strMsh.nodes(element, :);
    
    %% 2iii. Get the function handle for this type of loading
    if ischar(propNBC.fctHandle(iElmnt, :))
        loadFctHandle = str2func(propNBC.fctHandle(iElmnt, :));
    elseif isa(propNBC.fctHandle{iElmnt}, 'function_handle')
        loadFctHandle = propNBC.fctHandle{iElmnt};
    else
        error('Variable stored in NBC.fctHandle(%d,:) is neither a string nor it defines a function handle',iElmnt);
    end
    
    %% 2iv. Assemble the Element Freedome Table (EFT)
    EFT = zeros(numDOFsEl, 1);
    for iEFT = 1:numNodesEl
        EFT(iEFT) = element(iEFT);
    end
    
    %% 2v. Loop over all Gauss Points
    for iGP = 1:numGP
        %% 2v.1. Transform the Gauss Point from the parameter to the physical space
        xGP = (1 - GP(iGP))*x1/2 + (1 + GP(iGP))*x2/2;
        
        %% 2v.2. Compute the basis functions at the Gauss Point
        [basisFctOnGP, isInside] = computeCST2DBasisFunctions ...
            (nodes(1, 1:2), nodes(2, 1:2), nodes(3, 1:2), xGP(1, 1),xGP(1, 2));
        if ~isInside
            error('The computation of the natural coordinates in the CST triangle has failed');
        end
        
        %% 2v.3. Sort the basis functions into the basis functions array
        N = zeros(1, numDOFsEl);
        for iBF = 1:numNodesEl
            N(1, iBF) = basisFctOnGP(iBF, 1);
        end
        
        %% 2v.4. Compute the applied load onto the Gauss Point
        fluxOnGP = squeeze(loadFctHandle(xGP(1, 1), xGP(1, 2), xGP(1, 3), ...
            t, propNBC));
        
        %% 2v.5. Compute the determinant of the transformation from the physical space to the parameter space
        detJxxi = norm(x2 - x1)/2;
        
        %% 2v.6. Compute the element load vector on the Gauss point
        FElOnGP = N'*fluxOnGP*detJxxi*GW(iGP);
        
        %% 2v.7. Assemble the contribution to the global load vector via the EFT
        F(EFT) = F(EFT) + FElOnGP; 
    end
end

%% 3. Appendix
if strcmp(outMsg, 'outputEnabled')
    computationalTime = toc;
    fprintf('\nComutation of the load vector took %.2d seconds \n\n', computationalTime);
    fprintf('_________________Load Vector Computation Ended________________\n');
    fprintf('##############################################################\n\n\n');
end

end