function homDOFs = set_homogeneousDirichletBoundaryConditions...
    (homDOFs, patchCAD, propHomBC)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
% 
% Sets the homogeneous Dirichlet boundary conditions for a given patch.
%
%       Input :
%     homDOFs : Array with the preexisting numbering of the DOFs with
%               prescribed zero value
%    patchCAD : Cell array containing structures with fields,
%                   .CP : Control Point coordinates and weights
%   propHomBC : Structure containing information about the homogenerous
%               Dirichlet boundary conditions,
%                   .noCnd : Numnber of conditions
%             .xiExtension : Extensions of the BCs along xi
%            .etaExtension : Extensions of the BCs along eta
%              .directions : Cartesian directions along which the BCs are
%                            applied
%
%      Output :
%     homDOFs : Array with the updated numbering of the DOFs with
%               prescribed zero value
%
% Function layout :
%
% 0. Read input
%
% 1. Add the homogeneous Dirichlet boundary conditions
%
%% Function main body

%% 0. Read input
if ~isfield(propHomBC,'noCnd')
    error('BC must contain the field .noCnd which defines the number of conditions');
else
    if ~isnumeric(propHomBC.noCnd)
        error('BC.cnd must be of numeric type');
    else
        if floor(propHomBC.noCnd)~=propHomBC.noCnd
            error('BC.cnd must be of integer type');
        end
    end
end
if ~isfield(propHomBC,'xiExtension')
    error('BC must contain the field .xiExtension which defines the xi extensions of the BCs');
else
    if length(propHomBC.xiExtension) ~= propHomBC.noCnd
        error('Field BC.xiExtension must have the size of BC.noCnd');
    end
end
if ~isfield(propHomBC,'etaExtension')
    error('BC must contain the field .etaExtension which defines the eta extensions of the BCs');
else
    if length(propHomBC.etaExtension) ~= propHomBC.noCnd
        error('Field BC.etaExtension must have the size of BC.noCnd');
    end
end
if ~isfield(propHomBC,'directions')
    error('BC must contain the field .directions which defines the Cartesian directions along which to apply the BCs');
else
    if length(propHomBC.directions) ~= propHomBC.noCnd
        error('Field BC.etaExtension must have the size of BC.noCnd');
    else
        for i = 1:propHomBC.noCnd
            if length(propHomBC.directions{i}) > 3
                error('Field BC.directions must be an array of size 3');
            else
                for j = 1:length(propHomBC.directions{i})
                    if propHomBC.directions{i}(1,j) < 1 || propHomBC.directions{i}(1,j) > 3
                        error('Field BC.directions{i}(1,j) must have a value between 1 and 3');
                    end
                end
            end
        end
    end
end

%% 1. Add the homogeneous Dirichlet boundary conditions
for i = 1:propHomBC.noCnd
    xiExtension = propHomBC.xiExtension{i};
    etaExtension = propHomBC.etaExtension{i};
    directions = propHomBC.directions{i};
    for j = directions
        homDOFs = findDofs3D(homDOFs,xiExtension,etaExtension,directions(1,j),patchCAD.CP);
    end
end

end