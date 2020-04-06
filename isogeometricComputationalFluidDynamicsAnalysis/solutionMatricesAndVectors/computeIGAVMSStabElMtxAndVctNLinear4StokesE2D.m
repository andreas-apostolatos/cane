function [KEl, massMtxEl, FBodyEl] = ...
    computeIGAVMSStabElMtxAndVctNLinear4StokesE2D ...
    (dR, xiSpan, p, xi, Xi, etaSpan, q, eta, Eta, CP, parameters, ...
    computeBodyForces, tauM, tauC, Jxxi, Hxxi)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the stiffness and the mass matrix as well as the element load 
% vector due to body forces corresponding to the application of the 
% Variational Multiscale Stabilization method (VMS) to the linear 
% isogemetric Stokes problem in 2D.
% 
%           Input : 
%              dR : The vector containing the NURBS basis functions and up to 
%                   their second derivatives at the Gauss Point as follows:
%                   dR = 
%                   [R dR/dxi d^2R/dxi^2 dR/deta d^2R/deta/dxi d^2R/deta^2]
%  xiSpan,etaSpan : The knot span indices
%          xi,eta : The coordinates of the Gauss Points in the parameter
%                   space
%             p,q : polynomial degrees
%              CP : The set of Control Point Coordinates and weights
%      parameters : Flow parameters
% computeBodyFrc : The source force vector
%           tauM : Stabilization parameter for the momentum equilibrium
%           tauC : Stabilization parameter for the divergence-free 
%                  condition with respect to the velocity
%           Jxxi : Jacobian of the transformation from the physical space 
%                  the NURBS parametric space evaluated at the Gauss Point
%           Hxxi : The Hessian matrix at the Gauss Point
%    isTransient : Flag on whether transient analysis is encountered
%
%         Output :
%            KEl : The element stiffness matrix
%      massMtxEl : The element mass matrix
%        FBodyEl : The element body force vector
%
% Function layout :
%
% 0. Read input
%
% 1. Initialize all necessary matrices and vectors
%
% 2. Re-arrange the body force vector
%
% 3. Compute all necessary matrices and vectors by looping over all element contributions
% ->
%    3i. Transform the first derivatives of the basis functions from the NURBS parameter space to the physical space
%
%   3ii. Transform the second derivatives of the basis functions from the NURBS parameter space to the physical space
%
%  3iii. Compute the velocity basis functions matrix
%
%   3iv. Compute the pressure basis functions matrix
%
%    3v. Compute the derivatives of the basis functions related to the velocity gradient matrix
%
%   3vi. Compute the divergence of the basis functions matrix related to the velocity divergence matrix
%
%  3vii. Compute the Laplacian of the basis functions related to the velocities
%
% 3viii. Compute the matrix containing the gradient of the basis functions related to pressure
% <-
% 4. Compute the diffusion matrix at the Gauss Point
%
% 5. Compute the velocity-pressure coupling matrix at the Gauss Point
%
% 6. Compute the mass matrix
%
%    6i. Matrices related to diffusion nue*Laplacian(v)
%
%   6ii. Matrices related to pressure gradient Nabla(p)
%
%  6iii. Matrix related to the divergence-free condition (actual elimination of the saddle-point problem)
%
% 7. Compute the stabilization vectors of the right-hand side (stabilizing the right-hand side linear functional)
%
% 8. Compute the stabilization matrices of the left-hand side (stabilizing the left-hand side form)
%
% 9. Compute the element stabilized mass matrix
%
% 10. Compute the stabilized force vector due to source at the Gauss point
%
%% Function main body

%% 0. Read input

% The number of Control Points that affect the element under study
numNodesEl = (p + 1)*(q + 1);      

% The number of DoFs that affect the element under study
numDOFsEl = 3*numNodesEl;

% Initialize flag for the nature of the body forces
isBodyForceNumeric = 0;
isBodyForceIdenticallyZero = 0;
if isnumeric(computeBodyForces)
    isBodyForceNumeric = 1;
    if norm(computeBodyForces) == 0
        isBodyForceIdenticallyZero = 1;
    end
end

%% 1. Initialize all necessary matrices and vectors

% Initialize matrix containing the derivatives for each basis function
dRdx = zeros(numNodesEl, 2);

% Initialize tranformation matrix for the second
% derivatives
DJxxi = zeros(5, 5);

% Assign the values of DJxxi
% First block of DJxxi
for ki = 1:2
    for kj = 1:2
        DJxxi(ki, kj) = Jxxi(ki, kj)^2;
    end
end

% Complete first row
DJxxi(1, 3) = 2*Jxxi(1, 1)*Jxxi(1, 2);
DJxxi(1, 4) = Hxxi(1, 1)^2;
DJxxi(1, 5) = Hxxi(2, 1)^2;

% Complete second row
DJxxi(2, 3) = 2*Jxxi(2, 1)*Jxxi(2, 2);
DJxxi(2, 4) = Hxxi(1, 2)^2;
DJxxi(2, 5) = Hxxi(2, 2)^2;

% Complete third row
DJxxi(3, 1) = Jxxi(1, 1)*Jxxi(2, 1);
DJxxi(3, 2) = Jxxi(1, 2)*Jxxi(2, 2);
DJxxi(3, 3) = Jxxi(1, 1)*Jxxi(2, 2)+Jxxi(2, 1)*Jxxi(1, 2);
DJxxi(3, 4) = Hxxi(1, 3);
DJxxi(3, 5) = Hxxi(2, 3);

% Last block is the Jacobian matrix
for ki = 1:2
    for kj = 1:2
        DJxxi(ki + 3, kj + 3) = Jxxi(ki, kj);
    end
end

% Initialization of the matrix with the second derivatives
ddRddx = zeros(numNodesEl, 3);

% Initialize pressure basis functions matrix
RPressure = zeros(1, numDOFsEl);

% Initialize the derivatives of the basis functions related to the velocity 
% matrix
dudx = zeros(1, numDOFsEl);
dudy = zeros(1, numDOFsEl);
dvdx = zeros(1, numDOFsEl);
dvdy = zeros(1, numDOFsEl);

% Initialize the divergence of the basis functions matrix related to the
% velocities
dVelocities = zeros(1, numDOFsEl);

% Initialize the Laplacian of the basis functions related to the velocities
laplacianVelocities = zeros(3, numDOFsEl);

% Initialize the matrix containing the gradient of the basis functions
% related to pressure
gradientPressure = zeros(3, numDOFsEl);

% Initialize the velocities basis function matrix only if the source term
% (body forces) are non-zero
if isBodyForceNumeric || ~isBodyForceIdenticallyZero
    RVelocities = zeros(3, numDOFsEl);
end

%% 2. Re-arrange the body force vector
if isBodyForceNumeric && ~isBodyForceIdenticallyZero
    % Re-arrange the source vector to account also for the pressure
    % degrees of freedom
    bodyForcesRearranged = [computeBodyForces(1) computeBodyForces(2) 0]';
elseif ~isBodyForceNumeric
    % Find the Cartesian coordinates 
    P = computeCartesianCoordinatesOfAPointOnBSplineSurface ...
        (xiSpan, p, xi, Xi, etaSpan, q, eta, Eta, CP, dR(:, 1));
    x = P(1);
    y = P(2);

    % Evaluate the source vector at this location
    [bx, by] = computeBodyForces(x, y);
    bodyForcesRearranged = [bx by 0]';
end

%% 3. Compute all necessary matrices and vectors by looping over all element contributions
for i = 1:numNodesEl
    %% 3i. Transform the first derivatives of the basis functions from the NURBS parameter space to the physical space
    dRdx(i, :) = Jxxi\dR(i, [2 4])';
    
    %% 3ii. Transform the second derivatives of the basis functions from the NURBS parameter space to the physical space
    
    % Sort the basis functions and their derivatives in the parameter space
    %               dRActualSortedXiEta = 
    % [R d^2R/dxi^2 d^2R/deta^2 d^2R/dxi/deta dR/dxi dR/deta]
    dRActualSortedXiEta = [dR(i, [3 6 5]) dR(i, [2 4])];
    
    % Compute the basis functions and their derivatives in the physical
    % space
    dRActualSortedXY = DJxxi\dRActualSortedXiEta';
    
    % Isolate the second derivatives of the basis functions in the physical
    % space ddRddx = [d^2R/dX^2 d^2R/dY^2 d^2R/dX/dY]
    ddRddx(i,1:3) = dRActualSortedXY(1:3);
    
    %% 3iii. Compute the velocity basis functions matrix
    RVelocities(1,3*i - 2) = dR(i, 1);
    RVelocities(2,3*i - 1) = dR(i, 1);
    
    %% 3iv. Compute the pressure basis functions matrix
    RPressure(1, 3*i) = dR(i, 1);
    
    %% 3v. Compute the derivatives of the basis functions related to the velocity gradient matrix
    dudx(1,3*i - 2) = dRdx(i, 1);
    dudy(1,3*i - 2) = dRdx(i, 2);
    dvdx(1,3*i - 1) = dRdx(i, 1);
    dvdy(1,3*i - 1) = dRdx(i, 2);
    
    %% 3vi. Compute the divergence of the basis functions matrix related to the velocity divergence matrix
    dVelocities(1, 3*i - 2) = dRdx(i, 1);
    dVelocities(1,3*i - 1) = dRdx(i, 2);
    
    %% 3vii. Compute the Laplacian of the basis functions related to the velocities
    laplacianVelocities(1, 3*i - 2) = ddRddx(i, 1) + ddRddx(i, 2);
    laplacianVelocities(2, 3*i - 1) = ddRddx(i, 1) + ddRddx(i, 2);
    
    %% 3viii. Compute the matrix containing the gradient of the basis functions related to pressure
    gradientPressure(1, 3*i) = dRdx(i, 1);
    gradientPressure(2, 3*i) = dRdx(i, 2);
end

%% 4. Compute the diffusion matrix at the Gauss Point
De = parameters.nue*(dudx'*dudx + dudy'*dudy + dvdx'*dvdx + dvdy'*dvdy);

%% 5. Compute the velocity-pressure coupling matrix at the Gauss Point
Pe = dVelocities'*RPressure;

%% 6. Compute the mass matrix
massMtxEl = RVelocities'*RVelocities;

%% 6. Compute the stabilization matrices on the Gauss Point
%% 6i. Matrices related to diffusion nue*Laplacian(v)

% Diffusion-diffusion stabilization:
Dd = - parameters.nue^2*(laplacianVelocities'*laplacianVelocities);

% Diffusion-pressure stabilization:
Dp = parameters.nue*(laplacianVelocities'*gradientPressure);

% Diffusion-Mass stabilization
Dm = parameters.nue*laplacianVelocities'*RVelocities;

%% 6ii. Matrices related to pressure gradient Nabla(p)

% Pressure-diffusion stabilization: (Continuity equation)
Pd = - Dp';

% Pressure-pressure stabilization:
Pp = gradientPressure'*gradientPressure;

% Pressure-mass stabilization:
Pm = gradientPressure'*RVelocities;

%% 6iii. Matrix relateisTransientd to the divergence-free condition (actual elimination of the saddle-point problem)
Div = dVelocities'*dVelocities; % (Continuity equation)

%% 7. Compute the stabilization vectors of the right-hand side (stabilizing the right-hand side linear functional)
if isBodyForceNumeric || ~isBodyForceIdenticallyZero
    % Stabilization vector related to nue*Laplacian(v) (reaction):
    Fd = parameters.nue*laplacianVelocities';

    % Stabilization due to pressure gradient Nabla(p): (Continuity equation)
    Fp = gradientPressure';
end
       
%% 8. Compute the stabilization matrices of the left-hand side (stabilizing the left-hand side form)
KEl = De + Pe' - Pe + tauM*(Dd + Dp + Pd + Pp) + tauC*Div;

%% 9. Compute the element stabilized mass matrix
massMtxEl = massMtxEl + tauM*(Dm + Pm);

%% 10. Compute the stabilized force vector due to source at the Gauss point
if ~isBodyForceIdenticallyZero
    % Compute the force vector at the Gauss Point
    FBodyEl = (RVelocities' + tauM*(Fd + Fp))*bodyForcesRearranged;
else
    % Return an empty array
	FBodyEl = [];
end

end