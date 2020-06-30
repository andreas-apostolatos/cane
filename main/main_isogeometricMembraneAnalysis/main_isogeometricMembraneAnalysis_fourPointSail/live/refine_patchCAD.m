function patchCAD = refine_patchCAD...
    (patchCAD, num_xi, num_eta, dp, dq, outputEnabled)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
% 
% Refines the given NURBS patch using the ph method, that is, first the
% degree is elevated to the desirable order and then knots are inserted so
% that the continuity across the knots is kept as high as possible.
%
%           Input :
%        patchCAD : Unrefined patch containing the following fields,
%                 .p,.q : Polynomial orders
%              .Xi,.Eta : Knot vectors
%                   .CP : Control Point coordinates and weights
%              .isNURBS : Flag on whether the underlying basis is a NURBS
%                         or a B-Spline
%   outputEnabled : Allows printing of information in the command window if
%                   set as 'outputEnabled'
%
% num_xi, num_eta : Number of knots to add along the xi and eta directions
%          dp, dq : Incremental increase of the polynomial order
%
%          Output :
%        patchCAD : Refined patch containing the following fields,
%                 .p,.q : Polynomial orders
%              .Xi,.Eta : Knot vectors
%                   .CP : Control Point coordinates and weights
%              .isNURBS : Flag on whether the underlying basis is a NURBS
%                         or a B-Spline
%
% Function layout :
%
% 1. Perform degree elevation
%
% 2. Perform degree elevation
%
%% Function main body

%% 1. Perform degree elevation
[patchCAD.Xi, patchCAD.Eta, patchCAD.CP, patchCAD.p, patchCAD.q] =...
    degreeElevateBSplineSurface...
    (patchCAD.p, patchCAD.q, patchCAD.Xi, patchCAD.Eta, patchCAD.CP, dp, dq, outputEnabled);

%% 2. Perform degree elevation
[patchCAD.Xi,patchCAD.Eta,patchCAD.CP] = ...
    knotRefineUniformlyBSplineSurface...
    (patchCAD.p, patchCAD.Xi, patchCAD.q, patchCAD. Eta, patchCAD.CP, num_xi,num_eta, outputEnabled);

end