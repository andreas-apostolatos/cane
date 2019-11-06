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
function [tanStiffMtx,resVct] = computeRearrangedProblemMtrcs4MortarIGAMembrane...
    (constMtx,BSplinePatches,connections,tanStiffMtx,resVct)
%% Function main body
%
% Returns the re-arranged tangent stiffness matrix and residual vector of
% the problem when the mortar method is used for the multipatch coupling of
% membrane multipatches.
%
%                 Input :
%              constMtx : The set of the mortar tranformation matrices
%                         corresponding to each interface
%        BSplinePatches : Its an array of structures {patch1,patch2,...} 
%                         each of the patch structures containing the 
%                         following information
%                                 .p,.q: Polynomial degrees
%                              .Xi,.Eta: knot vectors
%                                   .CP: Control Points coordinates and 
%                                        weights
%                              .isNURBS: Flag on whether the basis is a 
%                                        NURBS or a B-Spline
%                             .homDOFs : The global numbering of the
%                                        DOFs where homogeneous Dirichlet
%                                        boundary conditions are applied
%                           .inhomDOFs : The global numbering of the
%                                        DOFs where homogeneous Dirichlet
%                                        boundary conditions are applied
%                     .valuesInhomDOFs : Prescribed values to the DOFs 
%                                        where homogeneous Dirichlet
%                                        boundary conditions are applied
%                          .domainDOFs : All DOFs of the patch but those
%                                        which are not located at any patch
%                                        interface
%                               FGamma : The boundary applied force vector
%                                        over the B-Spline patch
%                         .DOFNumbering: Numbering of the DOFs sorted into
%                                        a 3D array
%                           .parameters: material parameters of the shell
%                                  .int: On the numerical integration
%                                         .type : 'default' or 'user'
%                                        .xiNGP : No. of GPs along xi-
%                                                 direction for stiffness 
%                                                 entries
%                                       .etaNGP : No. of GPs along eta-
%                                                 direction for stiffness 
%                                                 entries
%                                 .xiNGPForLoad : No. of GPs along xi-
%                                                 direction for load 
%                                                 entries
%                                .etaNGPForLoad : No. of GPs along eta-
%                                                 direction for load 
%                                                 entries
%                                   .nGPForLoad : No. of GPs along boundary
%           connections : Define the connection between the patches:
%                            .No : Number of connections
%                     .xiEtaCoup : [patchID1 patchID2 xi12 eta12 xi21 eta21
%                                   ...      ...    ...   ...  ...   ...]
%           tanStiffMtx : The unmodified tangent stiffness matrix of the
%                         problem
%                resVct : The unmodified residual vector of the problem
%
% Function layout :
%
% 1. Loop over all the connections in the multipatch model and re-arrange the system by elimination of the Lagrange Multipliers DOFs
% ->
%    1i. Get the ID of the slave patch
%
%   1ii. Get the numbering of the necessary Lagrange Multilpiers, interface and domain DOFs
%
%  1iii. Get the mortar transformation matrix for the given patch connection
%
%   1iv. Re-arrange the stiffness entries by elimination of the Lagrange Multipliers DOFs
%
%    1v. Re-arrange the residual vector entries by elimination of the Lagrange Multipliers DOFs
% <-
%
%% Function main body

%% 1. Loop over all the connections in the multipatch model and re-arrange the system by elimination of the Lagrange Multipliers DOFs
for iConnections = 1:connections.No
    %% 1i. Get the ID of the slave patch
    idSlave = connections.mortar(iConnections,2);
    
    %% 1ii. Get the numbering of the necessary Lagrange Multilpiers, interface and domain DOFs
    
    % Master :
    % ________
    
    masterDOFs_C = connections.masterDOFs{iConnections};
    
    % Slave :
    % _______
    
    slaveDOFs_C = connections.slaveDOFs{iConnections};
%     slaveDOFs_Omega = BSplinePatches{idSlave}.domainDOFs;
    slaveDOFs_Omega = BSplinePatches{idSlave}.EFTPatches(~ismember(BSplinePatches{idSlave}.EFTPatches,slaveDOFs_C));
   
    %% 1iii. Get the mortar transformation matrix for the given patch connection
    TMortar = constMtx{iConnections};
    
    %% 1iv. Re-arrange the stiffness entries by elimination of the Lagrange Multipliers DOFs
    %
    %            | K_dd^(i)            K_dc^(i)                      0          |
    % K_mortar = | K_cd^(i) 2*K_cc^(i) + T'*K_cc^(i)*T T'*K_cd^(j) -T'*K_cc^(j) |
    %            |   0                 K_dc^(j)*T        K_dd^(j)    0          |
    %            |   0                -K_cc^(j)*T           0       K_cc^(j)    |
    %
    
    % Second row
    
    tanStiffMtx(masterDOFs_C,masterDOFs_C) = tanStiffMtx(masterDOFs_C,masterDOFs_C) ...
        + 2*TMortar'*tanStiffMtx(slaveDOFs_C,slaveDOFs_C)*TMortar;
    
    tanStiffMtx(masterDOFs_C,slaveDOFs_Omega) = tanStiffMtx(masterDOFs_C,slaveDOFs_Omega) ...
        + TMortar'*tanStiffMtx(slaveDOFs_C,slaveDOFs_Omega);
    
    tanStiffMtx(masterDOFs_C,slaveDOFs_C) = tanStiffMtx(masterDOFs_C,slaveDOFs_C) ...
        - TMortar'*tanStiffMtx(slaveDOFs_C,slaveDOFs_C);
    
    % Third row
    
    tanStiffMtx(slaveDOFs_Omega,masterDOFs_C) = tanStiffMtx(slaveDOFs_Omega,masterDOFs_C) ...
        + tanStiffMtx(slaveDOFs_Omega,slaveDOFs_C)*TMortar;
    
    tanStiffMtx(slaveDOFs_Omega,slaveDOFs_C) = ...
        zeros(length(slaveDOFs_Omega),length(slaveDOFs_C));
    
    % Fourth row

    tanStiffMtx(slaveDOFs_C,masterDOFs_C) = tanStiffMtx(slaveDOFs_C,masterDOFs_C) ...
        - tanStiffMtx(slaveDOFs_C,slaveDOFs_C)*TMortar;
    
    tanStiffMtx(slaveDOFs_C,slaveDOFs_Omega) = ...
        zeros(length(slaveDOFs_C),length(slaveDOFs_Omega));
    
    %% 1v. Re-arrange the residual vector entries by elimination of the Lagrange Multipliers DOFs
    %            |       r_d^(i)        |
    % r_mortar = | r_c^(i) + T'*r_c^(j) |
    %            |       r_d^(j)        |
    %            |          0           |    
    
    % Second row
    resVct(masterDOFs_C,1) = resVct(masterDOFs_C,1) + ...
        TMortar'*resVct(slaveDOFs_C,1);
    
    % Fourth row
%     resVct(slaveDOFs_C,1) = ...
%         tanStiffMtx(slaveDOFs_C,slaveDOFs_C)*(-TMortar*dHat(masterDOFs_C) + dHat(slaveDOFs_C));
	resVct(slaveDOFs_C,1) = zeros(length(slaveDOFs_C),1);
end

end