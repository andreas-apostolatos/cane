function [dHat, FComplete, rankD, condK, minEig, minElAreaSize] = ...
    solve_DDMPenaltyIGAKirchhoffLoveShellMultipatchesLinear ...
    (BSplinePatches, connections, propCoupling, solve_LinearSystem, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Gets the geometrical, technical, and coupling information for an 
% arbitrary number of Kirchhoff-Love shell patches together with the 
% constraint, loading  conditions and the connectivity information and 
% returns the displacement vector for each patch. The applied method is the 
% geometrically linear penalty method.
%
%                Input :
%       BSplinePatches : Structure containing all the information regarding 
%                        the connections between the multipatches
%          connections : Define the connection between the patches:
%                            .No : Number of connections
%                     .xiEtaCoup : [patchID1 patchID2 xi12 eta12 xi21 eta21
%                                 ...      ...    ...   ...  ...   ...]
%         propCoupling : Properties of the multipatch coupling
%                           .alphaD : penalty factor for the displacement
%                                     coupling
%                           .alphaR : penalty factor for the rotation
%                                     coupling
%                             .intC : On the integration of the coupling
%                                     interface
%   solve_LinearSystem : Function handle to the linear equation system
%                        solver
%               outMsg : Whether or not to output message on refinement 
%                       progress
%                       'outputEnabled' : enables output information
%   
%               Output :
%                 dHat : Array containing the displacement field of each 
%                        patch in the coupled system
%            FComplete : The complete force vector of the multipatch system
%                rankD : The rank deficiency of the linear system
%                condK : The condition number of the linear system
%               minEig : The minimum eigenvalue to the linear system
%        minElAreaSize : Array containing the minimum element area size for
%                        each patch
%
% Function Layout :
%
% 0. Read input
%
% 1. Get the number of DOFs for each patch and make an element freedom table for each Lagrange Multipliers field within the patch
%
% 2. Rearrange the vectors containing the global numbering of the DOFs where Dirichlet boundary conditions are applied to account for the coupled system
%
% 3. Create an element freedom table for each patch
%
% 4. Get a DOF numbering for each patch
%
% 5. Compute the constant part of the coupled stiffness matrix using the Penalty method
%
% 6. Solve the linear system
%
% 7. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('____________________________________________________________\n');
    fprintf('############################################################\n');
    fprintf('Statical Geometrically linear analysis for the domain\n');
    fprintf('decomposition of the isogeometric Kirchhoff-Love shell using \n');
    fprintf('the penalty algorithm has been initiated  \n');
    fprintf('\n');
    fprintf('Penalty factor for the displacements = %d\n', propCoupling.alphaD);
    fprintf('Penalty factor for the rotations = %d\n', propCoupling.alphaR);
    fprintf('____________________________________________________________\n');
    fprintf('\n');
    tic;
end

%% 0. Read input

% Define analysis type
analysis.type = 'isogeometricKirchhoffLoveShellAnalysis';

% Initialize dummy variables
dHatSaved = 'undefined';
dHatDotSaved = 'undefined';
dHatDDotSaved = 'undefined';
dHatDot = 'undefined';
dHatDDot = 'undefined';
propNLinearAnalysis = 'undefined';
computeUpdatedGeometry = 'undefined';
plot_IGANLinear = 'undefined';
massMtx = 'undefined';
dampMtx = 'undefined';
updateDirichletBCs = 'undefined';
propIDBC = 'undefined';
isReferenceUpdated = 'undefined';
isCosimulationWithEmpire = 'undefined';
propGraph = 'undefined';

% Steady-state analysis
propTransientAnalysis.computeProblemMtrcsTransient = 'undefined';
propTransientAnalysis.isStaticStep = true;

% The applied analysis is steady-state
t = 0;

% Tabulation
tab = '\t';

% Number of patches
numPatches = length(BSplinePatches);

% Total number of DOFs for the multi-patch system
numDOFs = 0;
for iPatches = 1:numPatches
    numDOFs = numDOFs + 3*BSplinePatches{iPatches}.noCPs;
end

% Initialize auxiliary variables
sizePrevious = 0;

% On the application of weak Dirichlet boundary conditions
isWeakDBC = zeros(numPatches, 1);
for iPatches = 1:numPatches
    if ~isempty(BSplinePatches{iPatches}.weakDBC)
        if isfield(BSplinePatches{iPatches}.weakDBC, 'noCnd')
            if BSplinePatches{iPatches}.weakDBC.noCnd > 0
                isWeakDBC(iPatches,1) = true;
            end
        end
    end
end

% Initialize the displacement vector
dHat = zeros(numDOFs,1);

%% 1. Get the number of DOFs for each patch and make an element freedom table for each Lagrange Multipliers field within the patch
for iPatches = 1:numPatches
    BSplinePatches{iPatches}.noDOFs = 3*BSplinePatches{iPatches}.noCPs;
    noDOFsPatchLM = 0;
    if isWeakDBC(iPatches, 1)
        if strcmp(BSplinePatches{iPatches}.weakDBC.method,'lagrangeMultipliers')
            for iCnd = 1:BSplinePatches{iPatches}.weakDBC.noCnd
                noDOFsPatchLMCnd = 3*length(BSplinePatches{iPatches}.weakDBC.lambda{iCnd}.CP(:,1));
                noDOFsPatchLM = noDOFsPatchLM + noDOFsPatchLMCnd;
                if iCnd == 1
                    index = 3*BSplinePatches{iPatches}.noCPs;
                else
                    index = BSplinePatches{iPatches}.weakDBC.lambda{iCnd-1}.EFT(length(BSplinePatches{iPatches}.weakDBC.lambda{iCnd-1}.EFT));
                end
                BSplinePatches{iPatches}.weakDBC.lambda{iCnd}.EFT = index + 1:index + noDOFsPatchLMCnd;
            end
        end
    end
    BSplinePatches{iPatches}.noDOFs = BSplinePatches{iPatches}.noDOFs + noDOFsPatchLM;
end

%% 2. Rearrange the vectors containing the global numbering of the DOFs where Dirichlet boundary conditions are applied to account for the coupled system
homDOFs = [];
inhomDOFs = [];
valuesInhomDOFs = [];
masterDOFs = [];
slaveDOFs = [];
for iPatches = 1:numPatches
    if iPatches ~= 1
        % Determine the number of the DOFs of the previous patches
        sizePrevious = sizePrevious + 3*BSplinePatches{iPatches - 1}.noCPs;
        
        % Add the numbering DOFs where homogeneous Dirichlet boundary
        % conditions are applied from the patch level
        homDOFsPatch = sizePrevious + BSplinePatches{iPatches}.homDOFs;
        homDOFs = mergesorted(homDOFs,homDOFsPatch);
        
        % Add the numbering DOFs where inhomogeneous Dirichlet boundary
        % conditions are applied from the patch level
        inhomDOFsPatch = sizePrevious + BSplinePatches{iPatches}.inhomDOFs;
        inhomDOFs = mergesorted(inhomDOFs,inhomDOFsPatch);
        
        % Add the prescribed values to the DOFs where inhomogeneous 
        % Dirichlet boundary conditions are applied from the patch level
        valuesInhomDOFsPatch = sizePrevious + BSplinePatches{iPatches}.valuesInhomDOFs;
        valuesInhomDOFs = mergesorted(valuesInhomDOFs,valuesInhomDOFsPatch);
        
        % Add the numbering of the master DOFs in the patch level
        masterDOFsPatch = sizePrevious + BSplinePatches{iPatches}.masterDOFs;
        masterDOFs = mergesorted(masterDOFs,masterDOFsPatch);
        
        % Add the numbering of the master DOFs in the patch level
        slaveDOFsPatch = sizePrevious + BSplinePatches{iPatches}.slaveDOFs;
        slaveDOFs = mergesorted(slaveDOFs, slaveDOFsPatch);
    else
        homDOFs = BSplinePatches{iPatches}.homDOFs;
        inhomDOFs = BSplinePatches{iPatches}.inhomDOFs;
        valuesInhomDOFs = BSplinePatches{iPatches}.valuesInhomDOFs;
        masterDOFs = BSplinePatches{iPatches}.masterDOFs;
        slaveDOFs = BSplinePatches{iPatches}.slaveDOFs;
    end
end

% Vector of the global numbering of the unconstrained DOFs
freeDOFs = 1:numDOFs;
freeDOFs(ismember(freeDOFs, homDOFs)) = [];

%% 3. Create an element freedom table for each patch
for iPatches = 1:numPatches
    if iPatches == 1
        BSplinePatches{iPatches}.EFTPatches = 1:3*BSplinePatches{iPatches}.noCPs;
    else
        BSplinePatches{iPatches}.EFTPatches = ...
            BSplinePatches{iPatches - 1}.EFTPatches(length(BSplinePatches{iPatches - 1}.EFTPatches)) + ...
            1:BSplinePatches{iPatches - 1}.EFTPatches(length(BSplinePatches{iPatches - 1}.EFTPatches)) + ...
            3*BSplinePatches{iPatches}.noCPs;
    end
end

%% 4. Get a DOF numbering for each patch
for iPatches = 1:numPatches
    % Get the number of Control Points in xi-direction
    nxi = length(BSplinePatches{iPatches}.CP(:, 1, 1));
    
    % Get the number of Control Points in eta-direction
    neta = length(BSplinePatches{iPatches}.CP(1, :, 1));
    
    % Initialize the DOF numbering array
    BSplinePatches{iPatches}.DOFNumbering = zeros(nxi, neta, 3);
    
    % Initialize counter
    k = 1;
    for cpj = 1:neta
        for cpi = 1:nxi
            BSplinePatches{iPatches}.DOFNumbering(cpi, cpj, 1) = k;
            BSplinePatches{iPatches}.DOFNumbering(cpi, cpj, 2) = k + 1;
            BSplinePatches{iPatches}.DOFNumbering(cpi, cpj, 3) = k + 2;

            % Update counter
            k = k + 3;
        end
    end
end

%% 5. Compute the constant part of the coupled stiffness matrix using the Penalty method
KConstantPenalty = computeConstantMtxForDDMPenaltyIGAThinStructure...
    (BSplinePatches, connections, numDOFs, propCoupling);

%% 6. Solve the linear system
[dHat, ~, ~, ~, FComplete, rankD, condK, minEig, ~, ~, minElAreaSize] = ...
    solve_IGALinearSystem ...
    (analysis, dHatSaved, dHatDotSaved, dHatDDotSaved, BSplinePatches, ...
    connections, dHat, dHatDot, dHatDDot, KConstantPenalty, massMtx, ...
    dampMtx, @computeStiffMtxAndLoadVctDDMPenaltyIGAKirchoffLoveShell, ...
    computeUpdatedGeometry, freeDOFs, homDOFs, inhomDOFs, ...
    valuesInhomDOFs, updateDirichletBCs, masterDOFs, slaveDOFs, ...
    solve_LinearSystem, t, propCoupling, propTransientAnalysis, ...
    propNLinearAnalysis, propIDBC, plot_IGANLinear, isReferenceUpdated, ...
    isCosimulationWithEmpire, tab, propGraph, outMsg);

%% 7. Appendix
if strcmp(outMsg,'outputEnabled')
    computationalTime = toc;
    fprintf('Static geometrically linear analysis took %.2d seconds \n\n', computationalTime);
    fprintf('___________________Linear Analysis Ended____________________\n');
    fprintf('############################################################\n\n\n');
    fprintf('\n');
end

end
