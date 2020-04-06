%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _______________________________________________________               %
%   _______________________________________________________               %
%                                                                         %
%   Technische Universität München                                        %
%   Lehrstuhl für Statik, Prof. Dr.-Ing. Kai-Uwe Bletzinger               %
%   _______________________________________________________               %
%   _______________________________________________________               %
%                                                                         %
%                                                                         %
%   Authors                                                               %
%   _______________________________________________________________       %
%                                                                         %
%   Dipl.-Math. Andreas Apostolatos    (andreas.apostolatos@tum.de)       %
%   Dr.-Ing. Roland Wüchner            (wuechner@tum.de)                  %
%   Prof. Dr.-Ing. Kai-Uwe Bletzinger  (kub@tum.de)                       %
%   _______________________________________________________________       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [KLinearEl,KNLinearEl,massMtxEl,FBodyEl] = computeIGAVMSStabElMtxAndVctNLinear4NSE2D...
    (dR,xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,parameters,bodyForces,tauM,tauC,upe,upVector,Jxxi,Hxxi)
%% Function documentation
%
% Returns the element stiffness, tangent stiffness and mass matrix as well
% as the element load vector due to body forces corresponding to the 
% application of the Variational Multiscale Stabilization method (VMS) to 
% the nonlinear isogemetric Navier-Stokes problem in 2D.
% 
%        Input : 
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
%      bodyForces : The source vector for the convection-diffusion-reaction
%                   equation
%            tauM : Stabilization parameter for the momentum equilibrium
%            tauC : Stabilization parameter for the divergence-free 
%                   condition with respect to the velocity
%             upe : The element discrete solution vector of the previous
%                   Newton iteration
%        upVector : The solution vector [u v p]' at the Gauss point from 
%                   the previous Newton iteration step
%            Jxxi : Jacobian of the transformation from the physical space the NURBS
%                   parametric space evaluated at the Gauss Point
%            Hxxi : The Hessian matrix at the Gauss Point
%
%          Output :
%       KLinearEl : The element linear stiffness matrix
%      KNLinearEl : The element nonlinear stiffness matrix
%       massMtxEl : The element mass matrix
%         FBodyEl : The element body force vector
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
%
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
%    3v. Compute the derivatives of the basis functions related to the velocity gradient matrix
%
%   3vi. Compute the divergence of the basis functions matrix related to the velocity divergence matrix
%
%  3vii. Compute the Laplacian of the basis functions related to the velocities
%
% 3viii. Compute the matrix containing the gradient of the basis functions related to pressure
%
%   3ix. Compute the B-operator matrix for the convection term
%
%    3x. Compute the B-operator matrices related to each row of the tangent convection matrix
%
%        3x.1. On the first B-operator matrix
%
%        3x.2. On the second B-operator matrix
%
%        3iv.3. Compute matrix due to the stabilization of the body forces-convection
%
% 4. Compute the diffusion matrix at the Gauss Point
%
% 5. Compute the convection matrix at the Gauss Point
%
% 6. Compute the convection-tangent matrix at the Gauss Point
%
% 7. Compute the velocity-pressure coupling matrix at the Gauss Point
%
% 8. Compute the stabilization matrices on the Gauss Point
%
%    8i. Matrices related to diffusion nue*Laplacian(v)
%
%   8ii. Matrices related to convection v.Gradient(v)
%
%  8iii. Matrices related to pressure gradient Nabla(p)
%
%   8iv. Matrix related to the divergence-free condition (actual elimination of the saddle-point problem)
%
% 9. Compute the mass matrix
%
% 10. Compute the stabilization matrices on the mass
%
%     10i. Mass-diffusion stabilization
%
%    10ii. Mass-convection stabilization
%
%   10iii. Mass-pressure stabilization
%
% 11. Compute the stabilization vectors of the right-hand side (stabilizing the right-hand side linear functional)
%
% 12. Compute the element stiffness matrix
%
% 13. Compute the element stabilized mass matrix
%
% 14. Compute the stabilized force vector due to source at the Gauss point
%
%% Function main body

%% 0. Read input

% The number of Control Points that affect the element under study
nNodesEl = (p+1)*(q+1);      

% The number of DoFs that affect the element under study
nDOFsEl = 3*nNodesEl;

% Initialize flag for the nature of the body forces
isBodyForceNumeric = 0;
isBodyForceIdenticallyZero = 0;
if isnumeric(bodyForces)
    isBodyForceNumeric = 1;
    if norm(bodyForces) == 0
        isBodyForceIdenticallyZero = 1;
    end
end

%% 1. Initialize all necessary matrices and vectors

% Initialize matrix containing the derivatives for each basis function
dRdx = zeros(nNodesEl,2);

% Initialize tranformation matrix for the second
% derivatives
DJxxi = zeros(5,5);

% Assign the values of DJxxi
% First block of DJxxi
for ki=1:2
    for kj=1:2
        DJxxi(ki,kj) = Jxxi(ki,kj)^2;
    end
end

% Complete first row
DJxxi(1,3) = 2*Jxxi(1,1)*Jxxi(1,2);
DJxxi(1,4) = Hxxi(1,1)^2;
DJxxi(1,5) = Hxxi(2,1)^2;

% Complete second row
DJxxi(2,3) = 2*Jxxi(2,1)*Jxxi(2,2);
DJxxi(2,4) = Hxxi(1,2)^2;
DJxxi(2,5) = Hxxi(2,2)^2;

% Complete third row
DJxxi(3,1) = Jxxi(1,1)*Jxxi(2,1);
DJxxi(3,2) = Jxxi(1,2)*Jxxi(2,2);
DJxxi(3,3) = Jxxi(1,1)*Jxxi(2,2)+Jxxi(2,1)*Jxxi(1,2);
DJxxi(3,4) = Hxxi(1,3);
DJxxi(3,5) = Hxxi(2,3);

% Last block is the Jacobian matrix
for ki=1:2
    for kj=1:2
        DJxxi(ki+3,kj+3) = Jxxi(ki,kj);
    end
end

% Initialization of the matrix with the second derivatives
ddRddx = zeros(nNodesEl,3);

% Initialize pressure basis functions matrix
RPressure = zeros(1,nDOFsEl);

% Initialize the derivatives of the basis functions related to the velocity 
% matrix
dudx = zeros(1,nDOFsEl);
dudy = zeros(1,nDOFsEl);
dvdx = zeros(1,nDOFsEl);
dvdy = zeros(1,nDOFsEl);

% Initialize the divergence of the basis functions matrix related to the
% velocities
dVelocities = zeros(1,nDOFsEl);

% Initialize the Laplacian of the basis functions related to the velocities
laplacianVelocities = zeros(3,nDOFsEl);

% Initialize the matrix containing the gradient of the basis functions
% related to pressure
gradientPressure = zeros(3,nDOFsEl);

% Initialize the B-operator related to the convection matrix
BoperatorConvection = zeros(3,nDOFsEl);

% Initialize the B-operator matrices related to each row of the tangent
% convection matrix
xiBoperatorTangentConvection = zeros(nDOFsEl,nDOFsEl);
etaBoperatorTangentConvection = zeros(nDOFsEl,nDOFsEl);

% Initialize the velocities basis function matrix only if the source term
% (body forces) are non-zero as well as the matrix due to the stabilization 
% of the body forces-convection term which goes on the left-hand side
if isBodyForceNumeric || ~isBodyForceIdenticallyZero
    RVelocities = zeros(3,nDOFsEl);
    Bc = zeros(nDOFsEl,nDOFsEl);
end

%% 2. Re-arrange the body force vector
if isBodyForceNumeric && ~isBodyForceIdenticallyZero
    % Re-arrange the source vector to account also for the pressure
    % degrees of freedom
    bodyForcesRearranged = [bodyForces(1) bodyForces(2) 0]';
elseif ~isBodyForceNumeric
    % Find the Cartesian coordinates 
    P = computeCartesianCoordinatesOfAPointOnBSplineSurface(xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,dR(:,1));
    x = P(1);
    y = P(2);

    % Evaluate the source vector at this location
    [bx,by] = bodyForcesForAnalyticalSolutionToStokesProblemInUnitSquare(x,y);
    bodyForcesRearranged = [bx by 0]';
end

%% 3. Compute all necessary matrices and vectors by looping over all element contributions
for i=1:nNodesEl
    %% 3i. Transform the first derivatives of the basis functions from the NURBS parameter space to the physical space
    dRdx(i,:) = Jxxi\dR(i,[2 4])';
    
    %% 3ii. Transform the second derivatives of the basis functions from the NURBS parameter space to the physical space
    
    % Sort the basis functions and their derivatives in the parameter space
    %               dRActualSortedXiEta = 
    % [R d^2R/dxi^2 d^2R/deta^2 d^2R/dxi/deta dR/dxi dR/deta]
    dRActualSortedXiEta = [dR(i,[3 6 5]) dR(i,[2 4])];
    
    % Compute the basis functions and their derivatives in the physical
    % space
    dRActualSortedXY = DJxxi\dRActualSortedXiEta';
    
    % Isolate the second derivatives of the basis functions in the physical
    % space ddRddx = [d^2R/dX^2 d^2R/dY^2 d^2R/dX/dY]
    ddRddx(i,1:3) = dRActualSortedXY(1:3);
    
    %% 3iii. Compute the velocity basis functions matrix
    RVelocities(1,3*i-2) = dR(i,1);
    RVelocities(2,3*i-1) = dR(i,1);   
    
    %% 3iv. Compute the pressure basis functions matrix
    RPressure(1,3*i) = dR(i,1);
    
    %% 3v. Compute the derivatives of the basis functions related to the velocity gradient matrix
    dudx(1,3*i-2) = dRdx(i,1);
    dudy(1,3*i-2) = dRdx(i,2);
    dvdx(1,3*i-1) = dRdx(i,1);
    dvdy(1,3*i-1) = dRdx(i,2);
    
    %% 3vi. Compute the divergence of the basis functions matrix related to the velocity divergence matrix
    dVelocities(1,3*i-2) = dRdx(i,1);
    dVelocities(1,3*i-1) = dRdx(i,2);
    
    %% 3vii. Compute the Laplacian of the basis functions related to the velocities
    laplacianVelocities(1,3*i-2) = ddRddx(i,1) + ddRddx(i,2);
    laplacianVelocities(2,3*i-1) = ddRddx(i,1) + ddRddx(i,2);
    
    %% 3viii. Compute the matrix containing the gradient of the basis functions related to pressure
    gradientPressure(1,3*i) = dRdx(i,1);
    gradientPressure(2,3*i) = dRdx(i,2);
    
    %% 3ix. Compute the B-operator matrix for the convection term
    BoperatorConvection(1,3*i-2) = upVector(1,1)*dRdx(i,1) + upVector(2,1)*dRdx(i,2);
    BoperatorConvection(2,3*i-1) = BoperatorConvection(1,3*i-2);
    
    %% 3x. Compute the B-operator matrices related to each row of the tangent convection matrix
    for j=1:nNodesEl
        %% 3x.1. On the first B-operator matrix
        xiBoperatorTangentConvection(3*j-2,3*i-2) = dR(j,1)*dRdx(i,1);
        xiBoperatorTangentConvection(3*j-1,3*i-2) = dR(j,1)*dRdx(i,2);

        %% 3x.2. On the second B-operator matrix
        etaBoperatorTangentConvection(3*j-2,3*i-1) = dR(j,1)*dRdx(i,1);
        etaBoperatorTangentConvection(3*j-1,3*i-1) = dR(j,1)*dRdx(i,2);     
        
        %% 3iv.3. Compute matrix due to the stabilization of the body forces-convection
        if isBodyForceNumeric && ~isBodyForceIdenticallyZero
            Bc(3*j-2,3*i-2) = - bodyForcesRearranged(1)*dR(i,1)*dRdx(j,1);
            Bc(3*j-2,3*i-1) = - bodyForcesRearranged(1)*dR(i,1)*dRdx(j,2);
            Bc(3*j-1,3*i-2) = - bodyForcesRearranged(2)*dR(i,1)*dRdx(j,1);
            Bc(3*j-1,3*i-1) = - bodyForcesRearranged(2)*dR(i,1)*dRdx(j,2);
        end
    end
end

%% 4. Compute the diffusion matrix at the Gauss Point
De = parameters.nue*(dudx'*dudx + dudy'*dudy + dvdx'*dvdx + dvdy'*dvdy);

% Auxiliary tangent convection B-operator matrix
CTangentBOperator = [upe'*xiBoperatorTangentConvection' 
                     upe'*etaBoperatorTangentConvection' 
                     zeros(1,nDOFsEl)];

%% 5. Compute the convection matrix at the Gauss Point
Ce = RVelocities'*BoperatorConvection;

%% 6. Compute the convection-tangent matrix at the Gauss Point
KNLinearEl = RVelocities'*CTangentBOperator;

%% 7. Compute the velocity-pressure coupling matrix at the Gauss Point
Pe = - dVelocities'*RPressure;

%% 8. Compute the stabilization matrices on the Gauss Point
%% 8i. Matrices related to diffusion nue*Laplacian(v)

% Diffusion-diffusion stabilization:
Dd = - parameters.nue^2*(laplacianVelocities'*laplacianVelocities);

% Diffusion-convection stabilization:
Dc = parameters.nue*laplacianVelocities'*BoperatorConvection;

% Diffusion-pressure stabilization:
Dp = parameters.nue*(laplacianVelocities'*gradientPressure);

%% 8ii. Matrices related to convection v.Gradient(v)

% Convection-diffusion stabilization:
Cd = - Dc';

% Convection-convection stabilization:
Cc = BoperatorConvection'*BoperatorConvection;

% Convection-pressure stabilization:
Cp = BoperatorConvection'*gradientPressure;

%% 8iii. Matrices related to pressure gradient Nabla(p)

% Pressure-diffusion stabilization: (Continuity equation)
Pd = - Dp';

% Pressure-convection stabilization: (Continuity equation)
Pc = Cp';

% Pressure-pressure stabilization: (Continuity equation)
Pp = gradientPressure'*gradientPressure;

%% 8iv. Matrix related to the divergence-free condition (actual elimination of the saddle-point problem)
Div = dVelocities'*dVelocities; % (Continuity equation)

%% 9. Compute the mass matrix
massMtxEl = RVelocities'*RVelocities;

%% 10. Compute the stabilization matrices on the mass
%% 10i. Mass-diffusion stabilization
Md = parameters.nue*laplacianVelocities'*RVelocities;

%% 10ii. Mass-convection stabilization
Mc = BoperatorConvection'*RVelocities;

%% 10iii. Mass-pressure stabilization
Mp = gradientPressure'*RVelocities;

%% 11. Compute the stabilization vectors of the right-hand side (stabilizing the right-hand side linear functional)
if isBodyForceNumeric && ~isBodyForceIdenticallyZero
    % Stabilization vector related to nue*Laplacian(v) (reaction):
    Fd = parameters.nue*laplacianVelocities';

    % Stabilization due to pressure gradient Nabla(p): (Continuity equation)
    Fp = gradientPressure';
end
       
%% 12. Compute the element stiffness matrix
if isBodyForceNumeric && ~isBodyForceIdenticallyZero
    % Linear stiffness matrix needed for the residual computation
    KLinearEl = De + Ce + Pe - Pe' + tauM*(Dd + Dc + Dp + Cd + Cc + Cp + Pd + Pc + Pp + Bc) + tauC*Div;
else
    % Linear stiffness matrix needed for the residual computation
    KLinearEl = De + Ce + Pe - Pe' + tauM*(Dd + Dc + Dp + Cd + Cc + Cp + Pd + Pc + Pp) + tauC*Div;
end

%% 13. Compute the element stabilized mass matrix
massMtxEl = massMtxEl + tauM*(Md + Mc + Mp);
    
%% 14. Compute the stabilized force vector due to source at the Gauss point
if isBodyForceNumeric && ~isBodyForceIdenticallyZero
    FBodyEl = (RVelocities' + tauM*(Fd + Fp))*bodyForcesRearranged;
else
	FBodyEl = [];
end

end