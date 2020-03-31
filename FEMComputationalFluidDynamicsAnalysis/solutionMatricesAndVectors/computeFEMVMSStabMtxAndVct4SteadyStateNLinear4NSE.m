function [K, resVct, minElEdgeSize] = ...
    computeFEMVMSStabMtxAndVct4SteadyStateNLinear4NSE...
    (propAnalysis, up, upSaved, upDot, upDotSaved, uMeshALE, ...
    DOFNumbering, fldMsh, F, loadFactor, propFldDynamics, t, ...
    parameters, computeBodyForces, gaussInt)
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
% corresponding to the nonlinear equation system occuring from the time
% discretization using the Bossak scheme and the space discretization using
% the classical finite element analysis of the Navier-Stokes equations.
%
% Reference :
%
% Ramon Codina, "A stabilized finite element method for generalized 
% stationary incompressible flows", Computer Methods in Applied Mechanics 
% and Engineering, Vol. 190 (2001), 2681-2706
%
% Implementation :
%
% KRATOS opensourse project, Riccardo Rossi
%
% Revision on the nonlinear equation system :
%                       
%                          G(u) * u = H                               (1)
%
% Nonlinear equation system (1) can be solved with the Newton method
% provided that the tangent of the system can be computed namely :
%
%               GTangent(u) = G(u) + dG(u)/du * u                     (2)
%
% In equation (2) the term dG(u)/du is obviously a third order tensor. By
% defining the residual of the equation system as :
%
%                   resVct(u) = G(u) * u - H                          (3)
%
% we iterate over all Newton iterations until convergence up to a given
% tolerance has been achieved, namely we solve iteratively the following 
% equation system:
%
%            GTangent(u^(i-1)) * du^i = - resVct(u^(i-1))            (4.1)
%                       u^i = u^(i-1) + du^i                         (4.2)
%
%             Input :
%      propAnalysis : Structure containing general properties of the
%                     analysis,
%                           .type : Analysis type
%                up : The discrete solution vector of the previous 
%                     nonlinear iteration step
%           upSaved : The discrete solution vector of the previous time 
%                     step
%             upDot : The rate of the discrete solution vector of the
%                     previous  nonlinear iteration step
%        upDotSaved : The rate of the discrete solution vector of the
%                     previous time step
%          uMeshALE : The mesh motion velocity field on the nodes of the
%                     mesh
%      DOFNumbering : The global numbering of the DOFs arranged in a
%                     3-dimentional array
%            fldMsh : Nodes and elements of the fluid mesh
%                 F : The boundary applied flux vector
%        loadFactor : Load factor in case more than one load steps are
%                     performed within each time step
%   propFldDynamics : Transient analysis parameters :
%                             .scheme : The time integration method
%                          .alphaBeta : (parameter for the Bossak scheme)
%                              .gamma : (parameter for the Bossak scheme)
%                                 .T0 : Start time of the simulation
%                               .TEnd : End time of the simulation
%                                 .nT : Number of time steps
%                                 .dt : Time step
%                 t : The current time of the transient simulation
%        parameters : Flow parameters
% computeBodyForces : Function handle to the computation of the body force
%                     vector
%                dt : The time step for the transient analysis. If the 
%                     analysis is statical choose dt = 0;
%               int : Structure responsible for the integration
%
%            Output :
%                 K : The linearized tangent stiffness matrix
%          KNLinear : The nonlinear stiffness matrix
%            resVct : The residual vector corresponding to the Newton 
%                     linearization of the Bossak discrete equation system 
%                     applied to the transient Navier-Stokes equations
%     minElEdgeSize : The minimum element edge size in the mesh
%
% Function layout :
%
% 0. Read input
%
% 1. Create the element freedom tables for all elements at once
%
% 2. Get the element discrete solution vector of the previous Newton iteration step
%
% 3. Get the dicrete mesh velocity vector (0 since we don't do ALE) GOT TO BE FIXED
%
% 4. Get the coordinates of the nodes in a matrix form
%
% 5. Get the minimum element edge size
%
% 6. Choose an integration rule
%
% 7. Loop over all the quadrature points
% ->
%    7i. Transform the Gauss Point location from the parameter to the physical space
%
%   7ii. Compute the basis functions and their derivatives at the Gauss Point
%
%  7iii. Compute the determinant of the Jacobian transformation from the physical to the parent space
%
%   7iv. Compute the tangent stiffness, the mass matrix and the body force vector on the Gauss Point
%
%    7v. Compute the load vector corresponding to the body force vector of the system
%
%   7vi. Add the contributions from the Gauss point
% <-
%
% 8. Add the contribution from the Gauss Point and assemble to the global system
%
% 9. Compute the system matrix corresponding to the Bossak time integration scheme
%
% 10. Compute the right-hand side vector corresponding to the Bossak time integration scheme
%
%% Function main body

%% 0. Read input

% Total number of nodes in the mesh
noNodes = length(fldMsh.nodes(:, 1));

% Number of DOFs per node
if strcmp(propAnalysis.type, 'NAVIER_STOKES_2D')
    noDOFsPerNode = 3;
    isAnalysis3D = false;
    noNodesEl = 3;
elseif strcmp(propAnalysis.type, 'NAVIER_STOKES_3D')
    noDOFsPerNode = 4;
    isAnalysis3D = true;
    noNodesEl = 4;
else
    error('Wrong analysis type specified')
end

% Total number of elements in the mesh
noElmnts = length(fldMsh.elements(:, 1));

% Total number of degrees of freedom
noDOFs = noDOFsPerNode*noNodes;

% Number of degrees of freedom per element
noDOFsEl = noDOFsPerNode*noNodesEl;

% Initialize arrays
KLineaEl = zeros(noElmnts, noDOFsEl, noDOFsEl);
KNLineaEl = zeros(noElmnts, noDOFsEl, noDOFsEl);
massMtxEl = zeros(noElmnts, noDOFsEl, noDOFsEl);
FBodyEl = zeros(noDOFsEl, 1);

% Compute a nessecary pre-factor for the Bossak time integration scheme
preFactor = (1 - propFldDynamics.alphaBeta)/propFldDynamics.gamma/ ...
    propFldDynamics.dt;

% Initialize the global body force vector
FBody = zeros(noDOFs, 1);

%% 1. Create the element freedom tables for all elements at once
EFT = zeros(noDOFsEl, noElmnts);
for iEFT = 1:noNodesEl
    for counterDOFsPerNode = 1:noDOFsPerNode - 1
        EFT(noDOFsPerNode*iEFT, :) = noDOFsPerNode*fldMsh.elements(:, iEFT)';
        EFT(noDOFsPerNode*iEFT - (noDOFsPerNode - counterDOFsPerNode), :) = ...
            EFT(noDOFsPerNode*iEFT, :) - (noDOFsPerNode - counterDOFsPerNode);
    end
end

%% 2. Get the element discrete solution vector of the previous Newton iteration step
upEl = up(EFT);

%% 3. Get the dicrete mesh velocity vector
if ~char(uMeshALE)
    uMeshALEEL = uMeshALE(EFT);
else
    uMeshALEEL = zeros(noDOFsEl, 1);
end

%% 4. Get the coordinates of the nodes in a matrix form

% define function to calculate euclidean norm
euclideanNorm = @(nodes) sqrt(nodes(:, 1, 1).^2 + nodes(:, 2, 1).^2 + nodes(:, 3, 1).^2);

% Minimum element edge size
if isAnalysis3D
    % Get the nodes of the mesh
    nodes1 = fldMsh.nodes(fldMsh.elements(:, 1), :);
    nodes2 = fldMsh.nodes(fldMsh.elements(:, 2), :);
    nodes3 = fldMsh.nodes(fldMsh.elements(:, 3), :);
    nodes4 = fldMsh.nodes(fldMsh.elements(:, 4), :);
    
    % get element sizes
    h = min( [ euclideanNorm(nodes1 - nodes2) euclideanNorm(nodes1 - nodes3) ...
               euclideanNorm(nodes2 - nodes3) euclideanNorm(nodes4 - nodes1) ...
               euclideanNorm(nodes4 - nodes2) euclideanNorm(nodes4 - nodes3)], ...
               [], 2);
else
    % Get the nodes of the mesh
    nodes1 = fldMsh.nodes(fldMsh.elements(:, 1), :);
    nodes2 = fldMsh.nodes(fldMsh.elements(:, 2), :);
    nodes3 = fldMsh.nodes(fldMsh.elements(:, 3), :);
    
	% get element sizes
    h = min( [ euclideanNorm(nodes1 - nodes2) euclideanNorm(nodes1 - nodes3) ...
               euclideanNorm(nodes2 - nodes3)], [], 2);
end

%% 5. Get the minimum element edge size
minElEdgeSize = min(h);

%% 6. Choose an integration rule

% Get the number of Gauss Points in xi and eta directions
if strcmp(gaussInt.type, 'default')
    noGP = 1;
elseif strcmp(gaussInt.type, 'user')
    noGP = gaussInt.domainNoGP;
end

% Get the Gauss Point coordinates a weights
if isAnalysis3D
    [GP, GW] = getGaussRuleOnCanonicalTetrahedron(noGP);
else
    [GP, GW] = getGaussRuleOnCanonicalTriangle(noGP);
end

%% 7. Loop over all the quadrature points
for iGP = 1:noGP
    %% 7i. Transform the Gauss Point location from the parameter to the physical space
    if isAnalysis3D
        xGP = GP(iGP, 1)*nodes1 + GP(iGP, 2)*nodes2 + GP(iGP, 3)*nodes3 + ...
        (1 - GP(iGP, 1) - GP(iGP, 2) - GP(iGP, 3))*nodes4;
    else
        xGP = GP(iGP, 1)*nodes1 + GP(iGP, 2)*nodes2 + (1 - GP(iGP, 1) - GP(iGP, 2))*nodes3;
    end

    %% 7ii. Compute the basis functions and their derivatives at the Gauss Point
    if isAnalysis3D
        [dN, area] = computeCST3DBasisFunctionsAndFirstDerivatives ...
            (nodes1, nodes2, nodes3, nodes4, xGP(:,1,:), xGP(:,2,:), xGP(:,3,:));
    else
        [dN, area] = computeCST2DBasisFunctionsAndFirstDerivatives ...
            (nodes1, nodes2, nodes3, xGP(:, 1, :), xGP(:, 2, :));
    end
    
    %% 7iii. Compute the determinant of the Jacobian transformation from the physical to the parent space
    detJxxi = 2.*area;
    
    %% 7iv. Compute the tangent stiffness, the mass matrix and the body force vector on the Gauss Point
    [KLineaElOnGP, KNLineaElOnGP, massMtxElOnGP, FBodyElOnGP] = ...
        computeFEMVMSStabElTangentStiffMtxMassMtxLoadVctNLinear4NSE ...
        (xGP(1, 1), xGP(1, 2), xGP(1, 3), t, upEl, uMeshALEEL, dN, ...
        computeBodyForces, parameters, h, propFldDynamics, isAnalysis3D);

    %% 7v. Compute the load vector corresponding to the body force vector of the system
    if norm(FBodyElOnGP) ~= 0
        FBody(EFT) = FBodyElOnGP(EFT) + FBodyElOnGP*GW(iGP)*detJxxi;
    end
    
    %% 7vi. Add the contributions from the Gauss point
    KLineaEl = KLineaEl + pstimes(KLineaElOnGP*GW(iGP), detJxxi);
    KNLineaEl = KNLineaEl + pstimes(KNLineaElOnGP*GW(iGP), detJxxi);
    if strcmp(propFldDynamics.timeDependence, 'TRANSIENT')
        massMtxEl = massMtxEl + pstimes(massMtxElOnGP*GW(iGP), detJxxi);
    end
    FBodyEl = FBodyEl + FBodyElOnGP;
end

%% 8. Add the contribution from the Gauss Point and assemble to the global system
if strcmp(propFldDynamics.timeDependence, 'TRANSIENT')
    [KLinear, KNLinear, massMtx] = assembleSparseMatricies ...
        (EFT, noDOFs, noDOFsEl, KLineaEl, KNLineaEl, massMtxEl);
elseif strcmp(propFldDynamics.timeDependence, 'STEADY_STATE')
    [KLinear, KNLinear] = assembleSparseMatricies ...
        (EFT, noDOFs, noDOFsEl, KLineaEl, KNLineaEl);
else
    error('wrong time dependence selected, see input file');
end

%% 9. Compute the system matrix corresponding to the Bossak time integration scheme
if strcmp(propFldDynamics.timeDependence, 'TRANSIENT')
    K = preFactor*massMtx + KLinear + KNLinear;
elseif strcmp(propFldDynamics.timeDependence, 'STEADY_STATE')
    K = KLinear + KNLinear;
end 

%% 10. Compute the right-hand side vector corresponding to the Bossak time integration scheme
if strcmp(propFldDynamics.timeDependence, 'TRANSIENT')
    resVct = (preFactor*massMtx + KLinear)*up - (FBody + loadFactor*F) - ...
            ((1-propFldDynamics.alphaBeta)/propFldDynamics.gamma/ ...
            propFldDynamics.dt)*massMtx*upSaved - ...
            ((1-propFldDynamics.alphaBeta)/propFldDynamics.gamma - 1)*...
            massMtx*upDotSaved;
elseif strcmp(propFldDynamics.timeDependence, 'STEADY_STATE')
    resVct = KLinear*up - (FBody + loadFactor*F);
else
    error('wrong time dependence selected, see input file');
end

end
