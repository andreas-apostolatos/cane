function Ke = computeElStiffMtxIGABernoulliBeamLinear(p,G,dG,dR,parameters)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation 
%
% computation of the Bernoulli element stiffness matrix evaluated at a 
% given Gauss Point.
%
%      Input :
%          p : The polynomial degree of the curve
%          G : The mid-line not normalized tangent vector
%         dG : The derivative of the mid-line tangent vector
%         dR : The matrix containing all the NURBS basis functions and
%              their first and second derivatives at the Gauss point
% parameters : The technical and geometrical parameters of the beam
%
%     Output :
%         Ke : The Bernoullli element stiffness matrix at the given
%              parametric location
%
% Function layout :
%
% 0. Read input
%
% 1. Compute normal strain epsilon
%
% 2. Compute change in the curvature kappa
%
% 3. Compute the geometrical and technical prefactors
%
% 4. Compute the virtual work contributions
%
% 5. Compute the element stiffness matrix on the Gauss Point
%
%% Function main body

%% 0. Read input

% length of the normal vector
L = norm(G);

% Initializations
epsilon = zeros(2*(p+1),1);
kappa = zeros(2*(p+1),1);

%% 1. Compute normal strain epsilon

% Loop over all the knot span contributions
for j = 1:p+1
    epsilon(1+(j-1)*2,1) = dR(j,2)*G(1,1);
    epsilon(2+(j-1)*2,1) = dR(j,2)*G(2,1);
end

%% 2. Compute change in the curvature kappa

% Necessary variables for the curvature
d1 = ( 1 / L^2 ) * G(1,1) * G(2,1);
d2 = ( 1 / L^2 ) * G(2,1)^2 - 1;
d3 = 1 - ( 1 / L^2 ) * G(1,1)^2;
alpha1 = ( dG(1,1) * d1 + dG(2,1) * d3 )/L;
alpha2 = ( dG(1,1) * d2 - dG(2,1) * d1 )/L;

% Loop over all the knot span contributions
for j = 1:p+1
   kappa(1+(j-1)*2,1) = - alpha1 * dR(j,2) + G(2,1) * dR(j,3)/L;
   kappa(2+(j-1)*2,1) = - alpha2 * dR(j,2) - G(1,1) * dR(j,3)/L;
end

%% 3. Compute the geometrical and technical prefactors

% Axial contribution
EAxial = parameters.EYoung*parameters.A / L^4;

% Bending contribution
EBending = parameters.EYoung*parameters.I / L^4;

%% 4. Compute the virtual work contributions

% virtual work done by the normal force n 
workEpsilon = EAxial * (epsilon*epsilon');

% virtual work done by the bending moment m 
workKappa = EBending *(kappa*kappa');

%% 5. Compute the element stiffness matrix on the Gauss Point
Ke = workEpsilon + workKappa;

end
