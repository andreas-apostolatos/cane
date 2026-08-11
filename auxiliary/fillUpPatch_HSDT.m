function BSplinePatch = fillUpPatch_HSDT...
    (analysis, p, Xi, q, Eta, CP, isNURBS, parameters, homDOFs, inhomDOFs,...
    valuesInhomDOFs, weakDBC, cables, NBC, varargin)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    [Your Name] (based on Andreas Apostolatos)
%
%% Function documentation
%
% Returns a B-Spline patch structure for HSDT analysis with 5 DOF per 
% control point: [u, v, w, θx, θy]
%
%                Input :
%             analysis : Structure containing the analysis type
%                  p,q : Polynomial degrees
%               Xi,Eta : Knot vectors
%                   CP : Control Point coordinates and weights
%              isNURBS : Whether the surface is NURBS or B-Spline
%           parameters : Material and geometric parameters (including 
%                        shear correction factor for HSDT)
%              homDOFs : Homogeneous Dirichlet boundary conditions
%            inhomDOFs : Inhomogeneous Dirichlet boundary conditions  
%     valuesInhomDOFs : Values for inhomogeneous Dirichlet BCs
%              weakDBC : Weak Dirichlet boundary conditions
%               cables : Embedded cable structures
%                  NBC : Neumann boundary conditions
%             varargin : Additional optional parameters
%
%               Output :
%         BSplinePatch : Complete B-Spline patch structure for HSDT
%
% Function layout :
%
% 1. Read input
%
% 2. Initialize the B-Spline patch structure
%
% 3. Fill up the patch with HSDT-specific information
%
%% Function main body

%% 1. Read input

% Number of control points
nxi = length(CP(:, 1, 1));
neta = length(CP(1, :, 1));
noCPs = nxi * neta;

% Handle optional input arguments
if nargin > 14
    int = varargin{5};
else
    % Default integration settings
    int.type = 'default';
end

%% 2. Initialize the B-Spline patch structure

% Basic patch information
BSplinePatch.p = p;
BSplinePatch.q = q;
BSplinePatch.Xi = Xi;
BSplinePatch.Eta = Eta;
BSplinePatch.CP = CP;
BSplinePatch.isNURBS = isNURBS;
BSplinePatch.noCPs = noCPs;

% Material parameters (must include shear correction factor)
BSplinePatch.parameters = parameters;
if ~isfield(parameters, 'shearCorrection')
    % Default shear correction factor for HSDT
    BSplinePatch.parameters.shearCorrection = 5/6;
    warning('Shear correction factor not specified. Using default value 5/6');
end

% Integration settings
BSplinePatch.int = int;

%% 3. Fill up the patch with HSDT-specific information

% Analysis type
BSplinePatch.analysis = analysis;

% Boundary conditions
BSplinePatch.homDOFs = homDOFs;
BSplinePatch.inhomDOFs = inhomDOFs;
BSplinePatch.valuesInhomDOFs = valuesInhomDOFs;
BSplinePatch.weakDBC = weakDBC;

% Embedded cables
if isempty(cables) || cables.No == 0
    BSplinePatch.cables.No = 0;
else
    BSplinePatch.cables = cables;
end

% Neumann boundary conditions
BSplinePatch.NBC = NBC;

% DOF numbering (5 DOF per control point)
% Will be set up in the solver function
BSplinePatch.DOFNumbering = [];

% Element freedom table
% Will be set up in the solver function  
BSplinePatch.EFTPatches = [];

% Total number of DOFs (5 per control point)
BSplinePatch.noDOFs = 5 * noCPs;

% HSDT-specific flags
BSplinePatch.isHSDT = true;
BSplinePatch.DOFsPerCP = 5;
BSplinePatch.DOFNames = {'u', 'v', 'w', 'theta_x', 'theta_y'};

% Initialize load vector (if needed for visualization)
if isfield(NBC, 'noCnd') && NBC.noCnd > 0
    BSplinePatch.FGamma = zeros(5 * noCPs, 1);
else
    BSplinePatch.FGamma = [];
end

% Additional HSDT information
BSplinePatch.shellTheory = 'HSDT';
BSplinePatch.description = 'Higher-order Shear Deformation Theory shell with 5 DOF per control point';

% Validation
if length(homDOFs) > 5 * noCPs
    error('Number of constrained DOFs exceeds total DOFs for HSDT formulation');
end

fprintf('HSDT patch initialized successfully:\n');
fprintf('  Control points: %d x %d = %d\n', nxi, neta, noCPs);
fprintf('  Total DOFs: %d (5 per control point)\n', 5 * noCPs);
fprintf('  Constrained DOFs: %d\n', length(homDOFs));
fprintf('  Free DOFs: %d\n', 5 * noCPs - length(homDOFs));
fprintf('  Shear correction factor: %.3f\n', BSplinePatch.parameters.shearCorrection);

end