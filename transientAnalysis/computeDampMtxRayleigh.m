function damMtxRayleigh = computeDampMtxRayleigh(constMtx,...
    computeProblemMatricesSteadyState,massMtx,noDOFs,BSplinePatches,...
    connections,propCoupling,t,propTransientAnalysis,noWeakDBCCnd,...
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
% Returns the damping matrix corresponding to the Rayleigh method.
%
%                               Input :
%                            constMtx : Constant to the problem matrix
%   computeProblemMatricesSteadyState : Function handle to the computation
%                                       of the problem matrix corresponding
%                                       to the steady-state problem
%                             massMtx : The mass matrix of the system
%                              noDOFs : Number of DOFs for the structure
%                      BSplinePatches : Array of B-Spline patches
%                                   .p,q : The polynomial orders of the 
%                                          B-Spline surface in both 
%                                          parametric directions
%                                .Xi,Eta : The knot vectors in both 
%                                          parametric directions
%                                    .CP : The set of control points and 
%                                          weights
%                               .isNURBS : Flag on whether the basis is a 
%                                          NURBS or a B-Spline
%                            .parameters : Technical parameters for the 
%                                          structure
%                               .homDOFs : The global numbering of the DOFs 
%                                          where homogeneous Dirichlet 
%                                          boundary conditions are applied
%                             .inhomDOFs : The global numbering of the DOFs 
%                                          where inhomogeneous Dirichlet 
%                                          boundary conditions are applied
%                       .valuesInhomDOFs : The values on the DOFs 
%                                          corresponding to the application 
%                                          of inhomogeneous Dirichlet 
%                                          boundary conditions
%                            masterDOFs : Master DOFs of the system
%                             slaveDOFs : Slave DOFs of the system
%                                  .NBC : Structure containing information 
%                                         on the application of the Neumann 
%                                         boundary conditions
%                                           .noCnd : Number of Neumann 
%                                                    boundary conditions
%                                 .xiLoadExtension : Cell array {.noCnd} 
%                                                    containing the load 
%                                                    extensions in the xi-
%                                                    direction
%                                .etaLoadExtension : Cell array {.noCnd} 
%                                                    containing the load 
%                                                    extensions in the eta-
%                                                    direction
%                                   .loadAmplitude : Array (1,.noCnd) 
%                                                    containing the load 
%                                                    amplitudes
%                                   .loadDirection : Array (1,.noCnd) 
%                                                    containing the load 
%                                                    directions
%                                  .computeLoadVct : Cell array {.noCnd} 
%                                                    containing the 
%                                                    function name for the 
%                                                    computation of the 
%                                                    load vector
%                                  .isConservative : Array (1,.noCnd) of 
%                                                    flags indicating 
%                                                    whether the load is 
%                                                    conservative or not               
%                    connections : Array defining the connection of the
%                                  B-Spline patches as well as the
%                                  extension of the coupling boundaries
%                   propCoupling : Coupling properties for a multipatch
%                                  geometry
%                              t : The time instance
%          propTransientAnalysis : On the transient analysis :
%                                   .damping.alpha : Alpha Rayleigh
%                                                    parameter
%                                    .damping.beta : Beta Rayleigh
%                                                    parameter
%                   noWeakDBCCnd : Number of weak Dirichlet boundary
%                                  conditions
%             isReferenceUpdated : Flag on whether the reference geometry 
%                                  is updated 
%                            tab : Tabulation for outputting information
%                                  onto the command window
%                         outMsg : Enables information outputting onto the
%                                  command window when chosen as 
%                                  'outputEnabled'
%
%                              Output :
%                      damMtxRayleigh : The damping matrix corresponding to
%                                       the Rayleigh approach
%
% Function layout :
%
% 0. Read input
%
% 1. Compute the linear stiffness matrix of the system
%
% 2. Compute the damping matrix corresppnding to the Rayleigh approach
%
%% Function main body
if ~isfield(propTransientAnalysis.damping,'alpha')
    error('Undefined Rayleigh damping parameter propTransientAnalysis.damping.alpha');
end
if ~isfield(propTransientAnalysis.damping,'beta')
    error('Undefined Rayleigh damping parameter propTransientAnalysis.damping.beta');
end
if ischar(massMtx)
    error('Undefined mass matrix');
end
if strcmp(outMsg,'outputEnabled')
    fprintf(strcat(tab,'Rayleigh parameter alpha = %f \n'),propTransientAnalysis.damping.alpha);
    fprintf(strcat(tab,'Rayleigh parameter beta = %f \n\n'),propTransientAnalysis.damping.beta);
end

%% 0. Read input

% Number of patches
noPatches = length(BSplinePatches);
if noPatches > 1
    noPatch = 'undefined';
else
    noPatch = 1;
end

% Define the dummy arrays
uSaved = 'undefined';
uDotSaved = 'undefined';
loadFactor = 'undefined';
noTimeStep = 'undefined';
propTransientAnalysisDummy = 'undefined';
if iscell(BSplinePatches)
    tanMtxLoad = struct([]);
    for iPatches = 1:length(BSplinePatches)
        tanMtxLoad{iPatches} = 'undefined';
    end
else
    tanMtxLoad = 'undefined';
end
noNonlinearIteration = 1;
u = zeros(noDOFs,1);
uDot = zeros(noDOFs,1);

%% 1. Compute the linear stiffness matrix of the system
[linStiffMtx,~,~,~,~] = ...
    computeProblemMatricesSteadyState...
    (constMtx,tanMtxLoad,u,uSaved,uDot,uDotSaved,BSplinePatches,...
    connections,propCoupling,loadFactor,noPatch,noTimeStep,...
    noNonlinearIteration,noWeakDBCCnd,t,propTransientAnalysisDummy,...
    isReferenceUpdated,tab,outMsg);

%% 2. Compute the damping matrix corresppnding to the Rayleigh approach
damMtxRayleigh = propTransientAnalysis.damping.alpha*massMtx + ...
    propTransientAnalysis.damping.beta*linStiffMtx;

end
