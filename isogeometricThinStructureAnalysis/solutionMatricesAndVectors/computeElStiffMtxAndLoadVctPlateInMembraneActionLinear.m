function Ke = computeElStiffMtxAndLoadVctPlateInMembraneActionLinear...
    (numCPsEl, dR, Jxxi, D)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Evaluates the stiffness matrix corresponding to the isogeometric plate in 
% membrane action of an element at the quadrature point
% 
%      Input :
%   noCPsLoc : Local number of Control Points
%         dR : The vector of the derivatives of the NURBS basis functions
%              at the quadrature point
%       Jxxi : jacobian of the transformation from the physical space the 
%              NURBS parametric space evaluated at the quadrature point
%          D : material matrix
%
%     Output :
%         Ke : The element stiffness matrix
%
% Function layout :
%
% 1. Compute the derivatives of the basis functions w.r.t. the physical domain at the quadrature point
%
% 2. Compute B-operator matrix at the parametric location where the basis functions are computed
%
% 3. Compute the element stiffness matrix at the parametric location where the basis functions are computed
%
%% Function main body

%% 1. Compute the derivatives of the basis functions w.r.t. the physical domain at the quadrature point

% Initialize matrix containing the derivatives for each basis function
dRdx = zeros(numCPsEl, 2);

% Compute the derivatives on the physical space given those at the parent
% domain
for i = 1:numCPsEl
    dRdx(i, :) = Jxxi\dR(i, :)';
end

%% 2. Compute B-operator matrix at the parametric location where the basis functions are computed

% Initialize B-operator matrix
B = zeros(3, 2*numCPsEl);

% Loop over all the entries
for i = 1:numCPsEl
   B(1,2*i-1) = dRdx(i,1);
   B(2,2*i) = dRdx(i,2);
   B(3,2*i-1) = dRdx(i,2);
   B(3,2*i) = dRdx(i,1);
end

%% 3. Compute the element stiffness matrix at the parametric location where the basis functions are computed
Ke = B'*D*B;

end
