function BSplinePatches = computeUpdatedGeometryIGAThinStructureMultipatches...
    (BSplinePatches, dHat)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Given the displacement field dHat corresponding to the solution of a
% multipatch system for the isogeometric Kirchhoff-Love shell, updates the
% displaced Control Point location.
%
%                    Input :
%           BSplinePatches : Its an array of structures {patch1,patch2,...}
%                            each of the patch structures containing the
%                            following information
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
%                               FGamma : The boundary applied force vector
%                                        over the B-Spline patch
%                           bodyForces : Function handle to the computation 
%                                        of the body forces
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
%                     dHat : The complete displacement field for the
%                            multipatch geometry
%
%                   Output :
%                            This function has no output arguments
%
% Function layout :
%
% 0. Read input
%
% 1. Update the displaced Control Point locations for each patch
% ->
%    1i. Get the vector of the unknown DOFs for the given patch
%
%   1ii. Get the actual Control Point displacement excluding Lagrange Multipliers DOFs
%
%  1iii. Update the Control Point coordinates of the patch
% <-
%
%% Function main body

%% 0. Read input

% Compute the number of multipatches in the system
noPatches = length(BSplinePatches);

%% 1. Update the displaced Control Point locations for each patch
for iPatches = 1:noPatches
    %% 1i. Get the vector of the unknown DOFs for the given patch
    dHatPatch = dHat(BSplinePatches{iPatches}.EFTPatches);
    
    %% 1ii. Get the actual Control Point displacement excluding Lagrange Multipliers DOFs
    nxi = length(BSplinePatches{iPatches}.CP(:,1,1));
    neta = length(BSplinePatches{iPatches}.CP(1,:,1));
    noDOFsDisp = 3*nxi*neta;
    dHatDisp = dHatPatch(1:noDOFsDisp);
    
    %% 1iii. Update the Control Point coordinates of the patch
    BSplinePatches{iPatches}.CPd = ...
        computeDisplacedControlPointsForIGAKirchhoffLoveShell...
        (BSplinePatches{iPatches}.CP,dHatDisp);
end

end
