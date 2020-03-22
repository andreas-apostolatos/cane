function [fldMsh] = correctionSolution(fldMsh,homDOFs,inhomDOFs,...
                    valuesInhomDOFs,propALE,solve_LinearSystem,...
                    propFldDynamics,i,u_flag)
%% Licensing
%
%  License:         BSD License
%                   cane Multiphysics default license: cane/license.txt
%
%  Main authors:    Matthew Keller
%
%% Function documentation
%
%  Adjusts the mesh depending on optimization output
%
%               Input :
%              fldMsh : Nodes and elements for the fluid mesh
%             homDOFs : The global numbering of the DOFs where homogeneous
%                       Dirichlet boundary conditions are applied
%           inhomDOFs : The global numbering of the DOFs where
%                       inhomogeneous Dirichlet boundary conditions are
%                       applied
%     valuesInhomDOFs : Prescribed values on the nodes where inhomogeneous
%                       Dirichlet boundary conditions are applied
%             propALE : Properties regarding the ALE boundary
%                         .nodes : The sequence of the nodal coordinates
%                                  on the ALE boundary
%                     .fcthandle : Function handle to the computation of
%                                  the ALE motion
%                       propUser : Extra user-defined parameters
%  solve_LinearSystem : Function handle to the solver for the linear 
%                       equation system
%     propFldDynamics : On the transient analysis :
%                             .method : The time integration method
%                          .alphaBeta : (parameter for the Bossak method)
%                              .gamma : (parameter for the Bossak method)
%                             .TStart : Start time of the simulation
%                               .TEnd : End time of the simulation
%                                 .nT : Number of time steps
%                                 .dt : Time step
%                   i : The current optimization iteration step
%              u_flag : Flag for ALE mesh motion
%
%              Output :
%              fldMsh : Updated mesh

%% Function main body

% Set perturbed conditions depending on user input
if u_flag == 1
    propALE.propUser.Perturb_Flag = 'dy'; % Mesh motion call to adjust height
elseif u_flag == 2
    propALE.propUser.Perturb_Flag = 'dx'; % Mesh motion call to adjust width
elseif u_flag == 3
    propALE.propUser.Perturb_Flag = 'dxdy'; % Mesh motion call to adjust taper ratio
elseif u_flag == 4
    propALE.propUser.Perturb_Flag = 'dtaper'; % Mesh motion call to adjust width
elseif u_flag == 5
    propALE.propUser.Perturb_Flag = 'dtaper_dy'; % Mesh motion call to adjust taper ratio
end

% Update the mesh for corrected state
[fldMsh,~,~,~] = computeUpdatedMeshAndVelocitiesPseudoStrALE2D...
    (fldMsh,homDOFs,inhomDOFs,valuesInhomDOFs,propALE,...
    solve_LinearSystem,propFldDynamics, i);

fldMsh.initialNodes = fldMsh.nodes;

%     graph.index = 1;  
%     graph.index = plot_referenceConfigurationFEMPlateInMembraneAction...
%     (fldMsh,'undefined',zeros(length(fldMsh.nodes(:,1)),1),[],graph,'');
%     graph.index = graph.index + 1;
end  