function plot_postprocResultantsIGAPlateInMembraneAction ...
    (p, q, Xi, Eta, CP, isNURBS, propParameters, xiGrid, etaGrid, ...
    dHat, propGraph)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Plots the reference configuration of an isogeometric plate in membrane 
% action together with the selected postprocessing resultant.
%
%          Input :
%            p,q : Polynomial degrees
%         Xi,Eta : Knot vectors in xi,eta-direction
%             CP : Control point coordinates and weights of the undeformed 
%                  plate
%        isNURBS : Flag on whether the geometrical basis is NURBS or 
%                  B-Spline
% propParameters : Structure containing information on the material
%                  parameters,
%                   .E : Young's modulus
%                 .nue : Poisson's ratio
%                    t : thickness
%           dHat : The displacement field of the control points
% xiGrid,etaGrid : The grid points used for the plotting of the NURBS
%                  geometry
%           dHat : The displacement field of the control points
%      propGraph : Structure containing information on the figures,
%                      .index : Index of the current figure
%
%         Output :
%                  graphics
%
% Function layout :
%
% 0. Read input
%
% 1. Compute the resultant array to be used for the visualization
%
% 2. Visualize the resultant and the knots over the domain
%
%% Function main body

%% 0. Read input

% Number of knots in xi-,eta-direction
numKnots_xi = length(Xi);
numKnots_eta = length(Eta);

% Assign a tolerance value
epsTol = 1e-9;

% Compute incremental steps for the resultant computation over the domain
% incremental step for eta:
deta = (Eta(numKnots_eta) - Eta(1))/(xiGrid - 1);

% incremental step for xi:
dxi = (Xi(numKnots_xi) - Xi(1))/(etaGrid-1);

% Initialize array of the resultant to be visualized
resultantComponent = zeros(xiGrid, etaGrid);

% Cartesian image of the parameter space
P = zeros(xiGrid, etaGrid, 3);

% Number of Control Points in xi-,eta-direction
numCPs_xi = length(CP(:, 1, 1));
numCPs_eta = length(CP(1, :, 1));

% Compute the material matrix for the plane stress problem
matMtx = propParameters.E/(1-propParameters.nue^2)*...
     [1                  propParameters.nue 0
     propParameters.nue 1                  0
     0                  0                  (1 - propParameters.nue)/2];

% Make a DOF numbering for the given patch
dHatElem = zeros(numKnots_xi - p - 1, numKnots_eta - q - 1, 2*(p + 1)*(q + 1));
for etaSpan = (q + 1):(numKnots_eta - q - 1)
    for xiSpan = (p+1):(numKnots_xi - p - 1)
        counter_xi = 1; 
        for c = etaSpan - q - 1:etaSpan - 1 
            for b = xiSpan - p:xiSpan
                dHatElem(xiSpan, etaSpan, counter_xi) = dHat(2*(c*numCPs_xi + b) - 1);
                dHatElem(xiSpan, etaSpan, counter_xi + 1) = dHat(2*(c*numCPs_xi + b));
                counter_xi = counter_xi + 2;
            end
        end
    end
end

% Initialize flags
isUndeformed = false;

%% 1. Compute the resultant array to be used for the visualization

% counting index in eta-direction
counter_eta = 1;  

% Initialize coordinate in eta-direction
eta = Eta(1);

% Loop over all the parametric coordinates in eta-direction
while eta <= Eta(numKnots_eta) + epsTol
    % Find the span in the eta-direction
    etaSpan = findKnotSpan(eta, Eta, numCPs_eta);
    
    % Initialize coordinate in xi-direction
    xi = Xi(1);
    
    % Initialize counter in xi-direction
    counter_xi = 1;
    
    % Loop over all the parametric coordinates in xi-direction
    while xi <= Xi(numKnots_xi) + epsTol
        % Find the span in xi-direction
        xiSpan = findKnotSpan(xi, Xi, numCPs_xi);
        
        % Compute the IGA basis functions and possibly their derivatives
        if strcmp(propGraph.resultant, 'displacement')
            numDrvs = 0;
        else
            numDrvs = 1;
        end
        dR = computeIGABasisFunctionsAndDerivativesForSurface ...
            (xiSpan, p, xi, Xi, etaSpan, q, eta, Eta, CP, isNURBS, numDrvs);
        
        % Compute the Cartesian image of the paratric point
        P(counter_xi, counter_eta, 1:3) = ...
            computeCartesianCoordinatesOfAPointOnBSplineSurface ...
            (xiSpan, p, xi, Xi, etaSpan, q, eta, Eta, CP, dR(:, 1));
        
        % Get the actual Control Point displacement vector for the current
        % knot span
        dHatActual(:, 1) = dHatElem(xiSpan, etaSpan, :);
        
        % Compute the resultant component
        if strcmp(propGraph.resultant, 'displacement')
            % Compute the displacement field at the given parametric
            % location
            d = computePostprocDisplacementIGAPlateInMembraneAction ...
                (p, q, dR(:, 1), dHatActual);
            
            % Decide upon the the visualization resultant
            if strcmp(propGraph.component, 'x')
                resultantComponent(counter_xi, counter_eta) = d(1, 1);
            elseif strcmp(propGraph.component,'y')
                resultantComponent(counter_xi, counter_eta) = d(2, 1);
            elseif strcmp(propGraph.component, '2norm')
                resultantComponent(counter_xi, counter_eta) = norm(d);
            end
        elseif strcmp(propGraph.resultant,'strain')
            % Compute the strain tensor in a Voigt notation
            epsilonVoigt = computePostprocVoigtStrainIGAPlateInMembraneAction...
                (xiSpan,p,etaSpan,q,CP,dR(:,2:3),dHatActual);
             
            % Compute the shear strain component, recall 
            % epsilon = [epsilonXX epsilonYY 2*epsilonXY]' 
            epsilonVoigt(3,1) = epsilonVoigt(3,1)/2;
            
            % Decide upon the the visualization resultant
            if strcmp(propGraph.component,'x')
                resultantComponent(counter_xi,counter_eta) = epsilonVoigt(1,1);
            elseif strcmp(propGraph.component,'y')
                resultantComponent(counter_xi,counter_eta) = epsilonVoigt(2,1);
            elseif strcmp(propGraph.component,'xy')
                resultantComponent(counter_xi,counter_eta) = epsilonVoigt(3,1);
            elseif strcmp(propGraph.component,'1Principal')
                resultantComponent(counter_xi,counter_eta) = 0.5*(epsilonVoigt(1,1)+epsilonVoigt(2,1) + sqrt((epsilonVoigt(1,1)-epsilonVoigt(2,1))^2 + epsilonVoigt(3,1)^2));
            elseif strcmp(propGraph.component,'2Principal')
                resultantComponent(counter_xi,counter_eta) = 0.5*(epsilonVoigt(1,1)+epsilonVoigt(2,1) - sqrt((epsilonVoigt(1,1)-epsilonVoigt(2,1))^2 + epsilonVoigt(3,1)^2));
            end
        elseif strcmp(propGraph.resultant,'stress')
            % Compute the strain tensor in a Voigt notation
            epsilonVoigt = computePostprocVoigtStrainIGAPlateInMembraneAction ...
                (xiSpan, p, etaSpan, q, CP, dR(:, 2:3), dHatActual);
            
            % Compute the stress tensor in a Voigt notation
            sigmaVoigt = matMtx*epsilonVoigt;
            
            % Decide upon the the visualization resultant
            if strcmp(propGraph.component, 'x')
                resultantComponent(counter_xi, counter_eta) = sigmaVoigt(1, 1);
            elseif strcmp(propGraph.component, 'y')
                resultantComponent(counter_xi, counter_eta) = sigmaVoigt(2, 1);
            elseif strcmp(propGraph.component, 'xy')
                resultantComponent(counter_xi, counter_eta) = sigmaVoigt(3, 1);
            elseif strcmp(propGraph.component, '1Principal')
                resultantComponent(counter_xi, counter_eta) = ...
                    0.5*(sigmaVoigt(1, 1) + sigmaVoigt(2,1) + sqrt((sigmaVoigt(1, 1) - sigmaVoigt(2, 1))^2 + 4*sigmaVoigt(3, 1)^2));
            elseif strcmp(propGraph.component, '2Principal')
                resultantComponent(counter_xi, counter_eta) = ...
                    0.5*(sigmaVoigt(1, 1) + sigmaVoigt(2, 1) - sqrt((sigmaVoigt(1, 1) - sigmaVoigt(2, 1))^2 + 4*sigmaVoigt(3, 1)^2));
            end
        end

        % Update the counter in xi-direction
        counter_xi = counter_xi + 1;
        
        % Update the parametric coordinate in xi-direction
        xi = xi + dxi;
    end
    % Update the counter in eta-direction
    counter_eta = counter_eta + 1;
    
    % Update the parametric coordinate in eta-direction
    eta = eta + deta;
end

%% 2. Visualize the resultant and the knots over the domain

% Plot the resultant over the reference configuration
surf(P(:, :, 1), P(:, :, 2), P(:, :, 3), resultantComponent(:, :));
hold on;

% Plot the element boundaries on the undeformed geometry
plot_knotsForBSplineSurfaceOnCartesianSpace ...
    (p, q, Xi, Eta, CP, isNURBS, isUndeformed, xiGrid, etaGrid);
hold off;

end