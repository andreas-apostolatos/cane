function Ke = computeElStiffMtxIGATimoshenkoBeamLinear(p,R,G,dG,parameters)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation 
%
% Returns the Timoshenko element stiffness matrix evaluated at the given
% parametric location for the case of the isogeometric Timoshenko beam
% analysis
%
%      Input :
%          p : The polynomial degree of the curve
%          R : The matrix containing all the NURBS basis functions and its
%              derivatives at the Gauss point
%          A : The mid-line not normalized tangent vector
%          H : The derivative of the mid-line tangent vector
% parameters : The technical and geometrical parameters of the beam
%
%     Output :
%         Ke : The element stiffness matrix for the isogeometric Timoshenko
%              beam analysis at the given parametric location
%
% Function layout :
%
% 0. Read input
%
% 1. Compute normal strain epsilon
%
% 2. Compute change in the curvature kappa
%
% 3. Compute shear strain gamma
%
% 4. Compute the geometrical and technical prefactors
%
% 5. Compute the virtual work contributions
%
% 6. Compute the element stiffness matrix on the Gauss Point
%
%% Function main body

% The length of the beam in the undeformed configuration
L = norm(G) ;

% Initializations
epsilon = zeros(3*(p+1),1);
kappa = zeros(3*(p+1),1);
gamma = zeros(3*(p+1),1);

%% 1. Compute normal strain epsilon

% Loop over all the knot span contributions
for j = 1:p+1
    epsilon(1+(j-1)*3,1) =  G(1,1) * R(j,2) ;
    epsilon(2+(j-1)*3,1) =  G(2,1) * R(j,2) ; 
    epsilon(3+(j-1)*3,1) = 0 ;
end

%% 2. Compute change in the curvature kappa
      
% Variable alpha
alpha = ( 1 / L^3 ) *( dG(1,1) * G(1,1) + dG(2,1) * G(2,1) ) ;

%Differentiation of A2 with respect to theta^1
A2xKomma1 = - ( 1 / L ) * dG(2,1) + alpha * G(2,1) ; 
A2yKomma1 = ( 1 / L ) * dG(1,1) - alpha * G(1,1) ;

% Variable b
beta = - G(1,1) * A2yKomma1 + G(2,1) * A2xKomma1 ;

%Computation of the curvature change kappa
for j = 1:p+1
    kappa(1+(j-1)*3,1) = A2xKomma1 * R(j,2) ; 
    kappa(2+(j-1)*3,1) = A2yKomma1 * R(j,2) ;
    kappa(3+(j-1)*3,1) = - L * R(j,2) + beta * R(j,2) ;
end

%% 3. Compute shear strain gamma
      
%Computation of the shear strain gamma
for j = 1:p+1
    gamma(1+(j-1)*3,1) = - G(2,1) * R(j,2) /  L ;
    gamma(2+(j-1)*3,1) = G(1,1) * R(j,2) / L ;
    gamma(3+(j-1)*3,1) = - L * R(j,1) ;
end

%% 4. Compute the geometrical and technical prefactors

% Axial contribution
EAxial = parameters.EYoung * parameters.A / L^4;

% Bending contribution
EBending = parameters.EYoung * parameters.I / L^4;

% Shear contribution
EShear = parameters.GShear * parameters.A / L^2;
      
%% 5. Compute the virtual work contributions

% Virtual work done by the normal force n 
epsilonWork = EAxial * (epsilon*epsilon');

% Virtual work done by the bending moment m 
kappaWork = EBending * (kappa*kappa');
      
% Virtual work done by the shear force q 
gammaWork = parameters.alpha * EShear * (gamma*gamma');
       
%% 6. Compute the element stiffness matrix on the Gauss Point
Ke = epsilonWork + kappaWork + gammaWork;

end
