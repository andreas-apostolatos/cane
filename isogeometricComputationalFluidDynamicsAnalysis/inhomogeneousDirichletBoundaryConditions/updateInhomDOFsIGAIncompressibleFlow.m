function [homDOFs, inhomDOFs, valuesInhomDOFs] = ...
    updateInhomDOFsIGAIncompressibleFlow ...
    (BSplinePatch, homDOFs, propIDBC, t)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the updated vector containing the global numbering of the DOFs
% where inhomogeneous Dirichlet boundary conditions are applied along with
% the vector of their prescribed time-dependent values. It is also removing
% corresponding DOF ids which are found in both the homogeneous and
% inhomogeneous Dirichlet boundary conditions.
%
%           Input :
%    BSplinePatch :
%        propIDBC : Structure containing information on the inhomogeneous 
%                   Dirichlet boundary conditions,
%               t : Time instance
%
%          Output :
%       inhomDOFs : The updated vector containing the global numbering of 
%                   DOFs where inhomogeneous Dirichlet boundary conditions 
%                   are applied
% valuesInhomDOFs : Vector containing the prescribed time-dependent values
%                   of the DOFs where inhomogeneous Dirichlet boundary 
%                   conditions are applied
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the inhomogeneous Dirichlet boundary conditions
% ->
%    1i. Get the Dirichlet boundary extension
%
%   1ii. Get the prescribed value
%
%  1iii. Get the direction of the flux
%
%   1iv. Determine whether this condition is unique over exisiting ones
%
%    1v. Compute the vector of the values on the DOFs where inhomogeneous Dirichlet boundary conditions are applied
%
% 2. Remove DOF ids which are found in both the homogeneous and the inhomogeneous Dirichlet boundary conditions
%
%% Function main body

%% 0. Read input

% Get the B-Spline patch
% Check the NURBS geometry input
if iscell(BSplinePatch)
    if length(BSplinePatch) > 1
        error('Multipatch NURBS surface is given as input to the computation of the stiffness matrix for a single patch NURBS surface');
    else
        BSplinePatch = BSplinePatch{1};
    end
end

% Get the Control Points of the patch
CP = BSplinePatch.CP;

% Get the flow parameters
parameters = BSplinePatch.parameters;

% Check if time dependent inhomogeneous Dirichlet boundary conditions are
% applied
if ~ischar(propIDBC)
    if ~isfield(propIDBC, 'numCnd')
        error('Field numCnd must be defined in structure propIDBC');
    end
    if ~isfield(propIDBC, 'prescribedValue')
        error('Field prescribedValue must be defined in structure propIDBC');
    end
    if ~isfield(propIDBC, 'prescribedDirection')
        error('Field prescribedDirection must be defined in structure propIDBC');
    end
    if ~isfield(propIDBC, 'isUniqueOnBoundary')
        error('Field isUniqueOnBoundary must be defined in structure propIDBC');
    end
    if ~isfield(propIDBC, 'prescribedValue')
        error('Field prescribedValue must be defined in structure propIDBC');
    else
        if isa(propIDBC.prescribedValue, 'function_handle')
            error('Field propIDBC.prescribedValue must define a function handle');
        end 
    end
end

% Intialize output arrays
inhomDOFs = [];
valuesInhomDOFs = [];
DOFsInBothHomAndInhomDBC = [];

%% 1. Loop over all the inhomogeneous Dirichlet boundary conditions and update their values
for iIDBC = 1:propIDBC.numCnd
    %% 1i. Get the Dirichlet boundary extension
    xib = propIDBC.xiExtension(iIDBC,:);
    etab = propIDBC.etaExtension(iIDBC,:);

    %% 1ii. Get the prescribed value
    prescribedValue = propIDBC.prescribedValue{iIDBC};

    %% 1iii. Get the direction of the flux
    prescribedDirection = propIDBC.prescribedDirection(iIDBC);

    %% 1iv. Determine whether this condition is unique over exisiting ones
    isUniqueOnBoundary = propIDBC.isUniqueOnBoundary(iIDBC);

    %% 1v. Compute the vector of the values on the DOFs where inhomogeneous Dirichlet boundary conditions are applied
    [inhomDOFs, valuesInhomDOFs] = ...
        computeVectorOfPreMotionOnInhomDoFsForIncompressibleFlow2D ...
        (inhomDOFs, valuesInhomDOFs, xib, etab, prescribedDirection, ...
        prescribedValue, isUniqueOnBoundary, CP, parameters, t);

end

%% 2. Remove DOF ids which are found in both the homogeneous and the inhomogeneous Dirichlet boundary conditions
if propIDBC.isDominant
    for iIDBC = 1:length(inhomDOFs)
        index = find(homDOFs == inhomDOFs(iIDBC));
        if norm(index) ~= 0
            DOFsInBothHomAndInhomDBC = ...
                mergesorted(DOFsInBothHomAndInhomDBC, index);
        end
    end
    homDOFs(DOFsInBothHomAndInhomDBC) = [];
else
    for iDBC = 1:length(homDOFs)
        index = find(inhomDOFs == homDOFs(iDBC));
        if norm(index) ~= 0
            DOFsInBothHomAndInhomDBC = ...
                mergesorted(DOFsInBothHomAndInhomDBC, index);
        end
    end
    inhomDOFs(DOFsInBothHomAndInhomDBC) = [];
end

end