function [u,uDot,uDDot,noTimeStep,propTransientAnalysis] = ...
    computeInitCndsFEMHeatTransferAnalysis...
    (analysis,strMsh,DOF4Output,parameters,propTransientAnalysis,...
    VTKResultFile,caseName,pathToOutput)
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
%
%           Input :
%        analysis : On the analysis,
%                   .type : Analysis type
%          strMsh : Nodes and elements in the mesh
%         homDOFs : The global numbering of the nodes where homogeneous
%                   Dirichlet boundary conditions are applied 
%       inhomDOFs : The global numbering of the nodes where inhomogeneous
%                   Dirichlet boundary conditions are applied
% valuesInhomDOFs : Prescribed values on the nodes where inhomogeneous
%                   Dirichlet boundary conditions are applied
%             NBC :    .nodes : The nodes where Neumann boundary 
%                                conditions are applied
%                    .loadType : The type of the load for each Neumann node
%                   .fctHandle : The function handle for each Neumann node
%                                for the computation of the load vector 
%                                (these functions are under the folder 
%                                load)
%  computeLoadVct : Function handle to the load vector computation
%               F : Global load vector corresponding to surface tractions
%      bodyForces : Function handle to body force vector computation
%        analysis : .type : The analysis type
%      parameters : Problem specific technical parameters
% nLinearAnalysis :     .scheme : The employed nonlinear scheme
%                    .tolerance : The residual tolerance
%                      .maxIter : The maximum number of the nonlinear 
%                                 iterations
%     strDynamics : .timeDependence : Steady-state or transient analysis
%                           .scheme : The time integration scheme
%                               .T0 : The start time of the simulation
%                             .TEnd : The end time of the simulation
%                          .nTSteps : The number of the time steps
%             int : On the numerical integration (quadrature)
%                      .type : 'default', 'manual'
%                       .nGP : Number of Gauss Points
%         intLoad : On the boundary integration
%                       .type : 'default', 'user'
%                       .nGP : No. of Gauss points%        
%        caseName : The name of the case in the inputGiD case folder
%    pathToOutput : Define the path to where to write out the results
%          outMsg : On outputting information
%
%          Output :
%            dHat : The nodal displacement field
%       minElSize : The minimum element area size in the mesh
%
% Function layout :
%
% 0. Read input
%
% 1. Find the prescribed and the free DOFs of the system
%
% 2. Solve the transient problem
%
% 3. Appendix
%
%% Function main body

% Get the number of nodes
numNodes = length(strMsh.nodes(:,1));

% Get the number of DOFs
numDOFs = numNodes;

% Assign the initial conditions
u = propTransientAnalysis.initialTemperature * ones(numDOFs,1);
uDot = zeros(numDOFs,1);
uDDot = zeros(numDOFs,1);


% Modify structure of transient analysis properties
% propTransientAnalysis.dt = propTransientAnalysis.noTimeSteps;
%propTransientAnalysis.dt = 0.1;
%propTransientAnalysis.dt =nTSteps
Tint = propTransientAnalysis.TEnd - propTransientAnalysis.T0;
propTransientAnalysis.dt = propTransientAnalysis.noTimeSteps;
propTransientAnalysis.noTimeSteps = Tint / propTransientAnalysis.dt;

noTimeStep = propTransientAnalysis.T0/propTransientAnalysis.dt;
% if(noTimeStep == 0)
%     noTimeStep = 1;
end
