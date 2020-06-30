function patchCAD = load_fourPointSinglePatchCAD(Length, Width, Height)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
% 
% Loads the basic CAD model of a four point geometry with alternating
% corner heights. This is the reference geometry of the four-point sail
% before form-finding.
%
%  Input :
% Length : Length of the geometry
%  Width : Width of the geometry
% Height : Alternating height of the four corner points
%
%      Output :
%    patchCAD : Cell array containing structures with fields,
%                 .p,.q : Polynomial orders
%              .Xi,.Eta : Knot vectors
%                   .CP : Control Point coordinates and weights
%              .isNURBS : Flag on whether the underlying basis is a NURBS
%                         or a B-Spline
%
% Function layout :
%
% 0. Read input
%
% 1. CAD modelling via NURBS
%
%% Function main body

%% 0. Read input
num_patches = 1;
patchCAD = cell(num_patches);

%% 1. CAD modelling via NURBS

% Polynomial degrees
patchCAD{1}.p = 1;
patchCAD{1}.q = 1;

% Knot vectors
patchCAD{1}.Xi = [0 0 1 1];
patchCAD{1}.Eta = [0 0 1 1];

% Control Point coordinates

% x-coordinates
patchCAD{1}.CP(:,:,1) = [-Length/2 Length/2
                         -Length/2 Length/2];
         
% y-coordinates
patchCAD{1}.CP(:,:,2) = [-Width/2 -Width/2
                         Width/2  Width/2];
         
% z-coordinates
patchCAD{1}.CP(:,:,3) = [Height 0    
                         0      Height];
       
% Weights
patchCAD{1}.CP(:,:,4) = [1 1
                         1 1];

% Find whether the geometrical basis is a NURBS or a B-Spline
patchCAD{1}.isNURBS = false;
nxi = length(patchCAD{1}.CP(:,1,1));
neta = length(patchCAD{1}.CP(1,:,1));
for i = 1:nxi
    for j = 1:neta
        if patchCAD{1}.CP(i,j,4)~=1
            patchCAD{1}.isNURBS = true;
            break;
        end
    end
    if patchCAD{1}.isNURBS
        break;
    end
end

end