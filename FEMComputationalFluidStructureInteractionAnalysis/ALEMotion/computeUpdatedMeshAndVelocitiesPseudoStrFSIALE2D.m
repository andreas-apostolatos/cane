function [fldMsh, uMeshALEFld, homDOFsFld, inhomDOFsFld, ...
    valuesInhomDOFsFld, freeDOFsFld] = ...
    computeUpdatedMeshAndVelocitiesPseudoStrFSIALE2D ...
    (fldMsh, nodesSaved, homDOFsFld, inhomDOFsFld, valuesInhomDOFsFld, ...
    freeDOFsFld, dHatInterface, propALE, propFSIFld, propFldDynamics)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the updated mesh node coordinates and the corresponding 
% velocities on the moving boundaries due the imposed motion on the FSI 
% interface. For the mesh motion the linear pseudo-structural
% solver is employed.
%
%              Input :
%             fldMsh : Structures containing the fluid mesh,
%                        .initialNodes : Initial nodes in the finite
%                                        element mesh
%                               .nodes : Current nodes in the finite
%                                        element mesh
%                            .elements : Elements in the finite element
%                                        mesh
%         nodesSaved : Array containing the IDs and coordinates of the 
%                      nodes from the previous time step
%         homDOFsFld : Global numbering of the DOFs on the fluid boundary
%                      where homogeneous Dirichlet boundary conditions are
%                      applied
%       inhomDOFsFld : Global numbering of the DOFs on the fluid boundary
%                      where inhomogeneous Dirichlet boundary conditions 
%                      are applied
% valuesInhomDOFsFld : Prescribed values of the DOFs on the fluid boundary
%                      where inhomogeneous Dirichlet boundary conditions 
%                      are applied
%        freeDOFsFld : Global numbering of the free DOFs regarding the
%                      fluid boundary value problem
%      dHatInterface : The prescribed mesh motion along the dry interface
%            propALE : Structure containing information about the ALE 
%                      motion of the fluid along the FSI interface,
%                      .nodes : IDs of the fluid nodes on the ALE boundary
%                  .fctHandle : Function handles to the computation of the 
%                               ALE motion at each node on the ALE boundary
%                     .isFree : Vector of flags indicating whether the 
%                               fluid motion is dictated by the ALE motion 
%                               at each node on the ALE boundary
%        propFSIFld : Structure containing information of the structural 
%                     wet surface,
%                       .nodes : Global numbering of the nodes on the wet 
%                                surface of the structure
%   propFldDynamics : Structure containing information on the fluid 
%                     dynamics,
%                                  .method : The time integration method
%                       .computeUpdatedVct : Function handle to the update
%                                            of the solution and its time
%                                            derivatives
%                               .alphaBeta : Parameter for the Bossak 
%                                            scheme
%                                   .gamma : Parameter for the Bossak 
%                                            scheme
%                                  .TStart : Start time of the simulation
%                                    .TEnd : End time of the simulation
%                             .noTimeSteps : Number of time steps
%                                      .dt : Time step size
%
%             Output :
%             fldMsh : The fluid mesh with the updated Cartesian
%                      coordinates of the nodes array,
%                       .nodes : Updated nodal coordinates due to the ALE
%                                motion
%        uMeshALEFld : Mesh velocity due to the ALE motion
%         homDOFsFld : Updated array due to the ALE motion containing the 
%                      global numbering of the DOFs where homogeneous 
%                      Dirichlet boundary conditions are applied
%       inhomDOFsFld : Updated array due to the ALE motion containing the 
%                      global numbering of the DOFs where inhomogeneous 
%                      Dirichlet boundary conditions are applied
% valuesInhomDOFsFld : Updated array due to the ALE motion containing the 
%                      prescribed values on the DOFs where inhomogeneous 
%                      Dirichlet boundary conditions are applied
%        freeDOFsFld : Updated array due to the ALE motion containing the
%                      free DOFs of the fluid boundary value problem
%
% Function layout :
%
% 0. Read input
%
% 1. Define the boundary conditions for the Arbitrary-Lagrangian Eulerian problem
%
% 2. Get the free DOFs of the Arbitrary-Lagrangian Eulerian problem
%
% 3. Solve the Arbitrary-Lagrangian Eulerian problem
%
% 4. Update the current Cartesian position of the nodes in the fluid mesh
%
% 5. Compute the velocity of the mesh regarding the Arbitrary-Lagrangian Eulerian problem
%
% 6. Get the global numbering of the DOFs and their prescribed values on the inhomogeneous Dirichlet boundary of the fluid problem on the dry surface
%
% 7. Update the global numbering of the DOFs and their proescribed values on the inhomogeneous Dirichlet boundary of the fluid problem 
%
% 8. Remove the global numbering of the DOFs where inhomogeneous Dirichlet boundary conditions are applied from the array containing the global numbering of the DOFs where homogeneous Dirichlet boundary conditions are applied
%
% 9. Update the array containing the global numbering of the free DOFs in the fluid problem
%
%% Function main body
if isempty(propALE) || ischar(propALE)
    error('propALE is undefined and mesh motion is called');
else
    if ~isfield(propALE, 'nodes')
        error('propALE must define a field with name nodes')
    end
    if ~isfield(propALE, 'fctHandle')
        error('propALE must define a field with name fctHandle')
    else
        for i = 1:length(propALE.fctHandle)
            fctHandle = propALE.fctHandle{i};
            if ~ischar(fctHandle)
                error('propALE.fctHandle{%d} must be a stringh defining a function handle', i);
            else
                if ~strcmp(fctHandle, 'computeNullH') && ...
                        ~strcmp(fctHandle, 'computeNullV') && ...
                        ~strcmp(fctHandle, 'computeNull0')
                    error('propALE.fctHandle{%d} = %s must be either "computeNullH", "computeNullV" or "computeNull0" for this function', ...
                        fctHandle);
                end
            end
        end
    end
    if ~isfield(propALE, 'isFree')
        error('propALE must define a field with name isFree')
    end
    if length(propALE.nodes(:, 1)) ~= length(propALE.fctHandle(:, 1)) || ...
            length(propALE.nodes(:, 1)) ~= length(propALE.isFree(:, 1)) || ...
            length(propALE.fctHandle(:, 1)) ~= length(propALE.isFree(:, 1))
        error('Arrays propALE.nodes, propALE.fctHandle and propALE.isFree must have the same length');
    end
end

%% 0. Read input

% Initialize dummy arrays
uSaved_pseudoStr = 'undefined'; 
uDotSaved_pseudoStr = 'undefined'; 
uDDotSaved_pseudoStr = 'undefined'; 
uDot_pseudoStr = 'undefined';
uDDot_pseudoStr = 'undefined';
massMtx_pseudoStr = 'undefined';
dampMtx_pseudoStr = 'undefined';
precompStiffMtx_pseudoStr = 'undefined';
precomResVct_pseudoStr = 'undefined';
propNLinearAnalysis_pseudoStr = 'undefined';
uMeshALEFld_pseudoStr = 'undefined';
propStrDynamics_pseudoStr = 'undefined';
t = 'undefined';
tab = 'undefined';
propInt_pseudoStr.type = 'default';
propAnalysis_pseudoStr.type = 'planeStress';

% Function handle to the solution of the linear equation syste
solve_LinearSystemALE = @solve_LinearSystemMatlabBackslashSolver;

% The material properties of the pseudo-structural solver
parameters_pseudoStr.nue = 0;
parameters_pseudoStr.E = 1e3;

% Get the number of nodes in the fluid computational domain
numNodes = length(fldMsh.nodes(:, 1));

% Initialize the vector of the mesh velocity due to the Arbitrary-
% Lagrangian Eulerian problem
uMeshALEFld = zeros(3*numNodes,1);

% Number of DOFs for the Arbitrary-Lagrangian Eulerian problem
numDOFsALE = 2*numNodes;

% Zero body forces
computeBodyForces_pseudoStr = @(x,y,z,t) zeros(length(x), 1, 3);

% Zero boundary forces
FGamma_pseudoStr = zeros(numDOFsALE, 1);

% Create a DOF numbering for the Arbitrary-Lagrangian Eulerian problem
DOFNumbering_presudoStr = 1:numDOFsALE;

% Initialize the solution of the Arbitrary-Lagrangian Eulerian problem
dHatALE = zeros(numDOFsALE, 1);

%% 1. Define the boundary conditions for the Arbitrary-Lagrangian Eulerian problem

% Get the global numbering of the DOFs from the boundary where homogeneous 
% Dirichlet boundary conditions are applied
[~, ~, nodeIdxALEFld] = intersect(propALE.nodes, fldMsh.nodes(:, 1), 'rows', 'stable');
homDOFs_pseudoStr = reshape(vertcat(2*nodeIdxALEFld' - 1, 2*nodeIdxALEFld'), 1, [])';

% Get the global numbering of the DOFs on the dry interface of the fluid
[~, ~, nodeIdxFld] = intersect(propFSIFld.nodes, fldMsh.nodes(:, 1), 'rows', 'stable');
inhomDOFs_preudoStr = reshape(vertcat(2*nodeIdxFld' - 1, 2*nodeIdxFld'), 1, [])';

% Get the values of the DOFs on the dry interface of the fluid
valuesInhomDOFs_pseudoStr = dHatInterface;

%% 2. Get the free DOFs of the Arbitrary-Lagrangian Eulerian problem
prescribedDOFs_pseudoStr = mergesorted(homDOFs_pseudoStr, inhomDOFs_preudoStr);
prescribedDOFs_pseudoStr = unique(prescribedDOFs_pseudoStr);
freeDOFs_pseudoStr = DOFNumbering_presudoStr;
freeDOFs_pseudoStr(ismember(freeDOFs_pseudoStr, prescribedDOFs_pseudoStr)) = [];

%% 3. Solve the Arbitrary-Lagrangian Eulerian problem
[dHatALE, ~, ~, ~] = solve_FEMLinearSystem...
    (propAnalysis_pseudoStr, uSaved_pseudoStr, uDotSaved_pseudoStr, ...
    uDDotSaved_pseudoStr, fldMsh, FGamma_pseudoStr, ...
    computeBodyForces_pseudoStr, parameters_pseudoStr, dHatALE, ...
    uDot_pseudoStr, uDDot_pseudoStr, massMtx_pseudoStr, ...
    dampMtx_pseudoStr, precompStiffMtx_pseudoStr, precomResVct_pseudoStr, ...
    @computeStiffMtxAndLoadVctFEMPlateInMembraneActionCST, ...
    DOFNumbering_presudoStr, freeDOFs_pseudoStr, ...
    homDOFs_pseudoStr, inhomDOFs_preudoStr, valuesInhomDOFs_pseudoStr, ...
    uMeshALEFld_pseudoStr, solve_LinearSystemALE, propStrDynamics_pseudoStr, ...
    t, propNLinearAnalysis_pseudoStr, propInt_pseudoStr, tab, '');

%% 4. Update the current Cartesian position of the nodes in the fluid mesh
fldMsh.nodes(:, 2) = fldMsh.initialNodes(:, 2) + dHatALE(1:2:end - 1, 1);
fldMsh.nodes(:, 3) = fldMsh.initialNodes(:, 3) + dHatALE(2:2:end, 1);

%% 5. Compute the velocity of the mesh regarding the Arbitrary-Lagrangian Eulerian problem
uMeshALEFld(1:3:end - 2, 1) = (fldMsh.nodes(:, 2) - nodesSaved(:, 2))/propFldDynamics.dt;
uMeshALEFld(1:3:end - 1, 1) = (fldMsh.nodes(:, 3) - nodesSaved(:, 3))/propFldDynamics.dt;

%% 6. Get the global numbering of the DOFs and their prescribed values on the inhomogeneous Dirichlet boundary of the fluid problem on the dry surface
inhomDOFsFld_FSI = reshape(vertcat(vertcat(3*nodeIdxFld' - 2, ...
    3*nodeIdxFld' - 1), 3*nodeIdxFld'), 1, [])';
valuesInhomDOFsFld_FSI = uMeshALEFld(inhomDOFsFld_FSI);

%% 7. Update the global numbering of the DOFs and their proescribed values on the inhomogeneous Dirichlet boundary of the fluid problem 
inhomDOFsFld = horzcat(inhomDOFsFld, inhomDOFsFld_FSI');
[inhomDOFsFld, idxSort] = sort(inhomDOFsFld);
valuesInhomDOFsFld = horzcat(valuesInhomDOFsFld, valuesInhomDOFsFld_FSI');
valuesInhomDOFsFld = valuesInhomDOFsFld(idxSort);

%% 8. Remove the global numbering of the DOFs where inhomogeneous Dirichlet boundary conditions are applied from the array containing the global numbering of the DOFs where homogeneous Dirichlet boundary conditions are applied
homDOFsFld(ismember(homDOFsFld, inhomDOFsFld)) = [];

%% 9. Update the array containing the global numbering of the free DOFs in the fluid problem
prescribedDOFsFld = mergesorted(homDOFsFld, inhomDOFsFld);
freeDOFsFld(ismember(freeDOFsFld, prescribedDOFsFld)) = [];

end