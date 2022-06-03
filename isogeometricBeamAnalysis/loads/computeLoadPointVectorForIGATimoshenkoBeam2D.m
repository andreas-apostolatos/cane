function [F, tanMtxLoad] = computeLoadPointVectorForIGATimoshenkoBeam2D ...
    (Fold, BSplinePatch, xib, etab, loadAmplitude, loadDir, isFollower, t, ...
    analysis, dir)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% returns the consistent nodal forces to a point load loadAmplitude for the 
% isogeometric 2D beams of type Bernoulli and Timoshenko.
%
%        Input :
%         Fold : The outdated force vector
% BSplinePatch : Structure containing all computational information for the
%                B-Spline patch,
%                   .p : Polynomial degree of the curve
%                  .Xi : Knot vector of the curve
%                  .CP : Control point coordinates and weights of the curve
%             .isNURBS : Flag on whether the underlying basis is a NURBS or
%                        a  B-Spline
%          xib : Load extension (e.g. [0.5 0.8])
%         etab : Dummy variable for this function
%      loadFnc : Load magnitude or function handle to its computation
%      loadDir : Load direction, 1=tangent 2=normal direction
%   isFollower : Dummy variable for this function
%            t : Time instance
%          int : Gaussian integration parameters
%       outMsg : Whether or not to output message on refinement progress
%                'outputEnabled' : enables output information
%
%        Output :
%             F : updated force vector
%    tanMtxLoad : Dummy output for this function
%
% 0. Read input
%
% 1. Evaluate the NURBS basis functions at the point of the load application
%
% 2. Compute the load component at the application point
%
% 3. If the given old load vector is not null, sum them up with the newly computed one
%
%% Function main body

%% 0. Read input

% Get the computational information from the patch
p = BSplinePatch.p;
Xi = BSplinePatch.Xi;
CP = BSplinePatch.CP;
isNURBS = BSplinePatch.isNURBS;

% Number of Control points
numCPs_xi = length(CP(:, 1));

% Find the knot span where the location of the load application is
knotSpan = findKnotSpan(xib, Xi, numCPs_xi);

% Initialize the load vector
F = zeros(size(CP,1)*3,1);

% Initialize dummy output
tanMtxLoad = 'undefined';

%% 1. Evaluate the NURBS basis functions at the point of the load application
R = computeIGABasisFunctionsAndDerivativesForCurve ...
    (knotSpan, p, xib, Xi, CP, isNURBS, 0);

%% 2. Compute the load component at the application point

% Initialize counter
k = 1;

% Loop over all the contributions at the current knot span
    
for b = 1:p + 1
    % Find the correct index for the current knot span
    index = dir + (knotSpan - p - 1)*2 + (b - 1)*2;
    
    % Compute the entries of the load vector according to the respective
    % theory
    index_tim = index + (knotSpan - p - 1) + b - 1;
    F(index_tim,1) = R(k)*loadAmplitude; 

    % Update counter
    k = k + 1;
end
    
%% 3. If the given old load vector is not null, sum them up with the newly computed one
if isvector(Fold)  
    F = F + Fold;   
end

end
