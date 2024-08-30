function epsilon = computePostprocVoigtStrainIGAPlateInMembraneAction...
    (xiSpan, p, etaSpan, q, CP, dR, dHatEl)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns strain vector epsilon = [epsilonXX epsilonYY 2*epsilonXY]' 
% corresponding to the isogeometric plate in membrane action at the given 
% parametric location.
%
%          Input : 
% xiSpan,etaSpan : The knot span indices
%            p,q : polynomial degrees in xi-,eta-directions
%             CP : Control Point coordinates and weights
%             dR : The first derivatives of the basis funnctions 
%         dHatEl : The Control Point displacement vector affecting the
%                  element where the basis functions are computed
%
%  Output :
% epsilon : element strain vector strain vector epsilon = [exx eyy 2exy]'
%
% Function layout :
%
% 0. Read input 
%
% 1. Compute the Jacobian of the transformation from the physical to the NURBS parameter space
%
% 2. Compute the derivatives of the basis functions w.r.t. the physical domain at the quadrature point
%
% 3. Compute B-operator matrix at the quadrature point
%
% 4. Compute the strain vector in a Voigt notation epsilon = [epsilonXX epsilonYY 2*epsilonXY]
%
%% Function main body

%% 0. Read input     

% Local number of Control Points
numCPsLoc = (p + 1)*(q + 1);

% Jacobian of transformation from the physical space to the NURBS parent
% domain
Jxxi = zeros(2, 2);

% Initialize counter
k = 0;

% Assign a tolerance value
tolerance = 1e-10;

%% 1. Compute the Jacobian of the transformation from the physical to the NURBS parameter space
for c = 0:q
    for b = 0:p
        % Update counter
        k = k + 1;
        Jxxi(1, 1) = Jxxi(1, 1) + CP(xiSpan - p + b, etaSpan - q + c, 1)*dR(k, 1);
        Jxxi(1, 2) = Jxxi(1, 2) + CP(xiSpan - p + b, etaSpan - q + c, 2)*dR(k, 1);
        Jxxi(2, 1) = Jxxi(2, 1) + CP(xiSpan - p + b, etaSpan - q + c, 1)*dR(k, 2);
        Jxxi(2, 2) = Jxxi(2, 2) + CP(xiSpan - p + b, etaSpan - q + c, 2)*dR(k, 2);
    end
end

%% 2. Compute the derivatives of the basis functions w.r.t. the physical domain at the quadrature point
if det(Jxxi) > tolerance
    % Initialize matrix containing the derivatives for each basis function
    dRdx = zeros(numCPsLoc, 2);

    % Compute the derivatives on the physical space given those at the parent
    % domain
    for i = 1:numCPsLoc
        dRdx(i, :) = Jxxi\dR(i, :)';
    end
end

%% 3. Compute B-operator matrix at the quadrature point
if det(Jxxi) > tolerance
    % Initialize B-operator matrix
    B = zeros(3, 2*numCPsLoc);

    % Loop over all the entries
    for i = 1:numCPsLoc
       B(1, 2*i - 1) = dRdx(i, 1);
       B(2, 2*i) = dRdx(i, 2);
       B(3, 2*i - 1) = dRdx(i, 2);
       B(3, 2*i) = dRdx(i, 1);
    end
end

%% 4. Compute the strain vector in a Voigt notation epsilon = [epsilonXX epsilonYY 2*epsilonXY]
if det(Jxxi) > tolerance
    epsilon = B*dHatEl;
else
    epsilon = zeros(3, 1);
end

end