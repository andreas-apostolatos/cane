function [epsilon, sigma] = ...
    computePostprocFEMPlateInMembraneActionCSTNLinear...
    (strMsh, analysis, parameters, dHat)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the strain field [epsilonXX epsilonYY epsilonXY] and stress field
% [sigmaXX sigmaYY sigmaXY] corresponding to a linear classical finite 
% element plate in membrane action analysis using the constant strain
% triangle.
%
%       Input :
%      strMsh : Nodes and elements in the mesh
%    analysis : .type : The analysis type
%  parameters : Problem specific technical parameters
%        dHat : The nodal displacement field
%
%      Output :
%     epsilon : The strain field [epsilonXX epsilonYY epsilonXY]
%       sigma : The stress field [sigmaXX sigmaYY sigmaXY]
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the elements in the mesh
% ->
%    1i. Get the element in the mesh
%
%   1ii. Cast the element with respect to the number of its nodes
%
%  1iii. Get the nodes in the element
%
%   1iv. Compute the basis functions and their derivatives at the Gauss Point
%
%    1v. Compute the displacement gradient components
%
%   1vi. Compute the strain vector in a Voigt notation at the current element
%
%  1vii. Compute the stress vector in a Voigt notation at the current element
%
% 1viii. Correct the shear component of the strain epsilon = [epsilonXX epsilonYY 2*epsilonXY]
% <-
%
%% Function main body

%% 0. Read input

% Number of nodes in the mesh
numElements = length(strMsh.elements);

% Compute the material matrix for the given problem (The shear entry is 
% multiplied by two so that it returns the true strain field and not the 
% one needed for the computation of the internal virtual work)
if strcmp(analysis.type, 'planeStress')
    preFactor = parameters.E/(1 - parameters.nue^2);
    C = preFactor*[1              parameters.nue  0
                   parameters.nue 1               0
                   0              0               (1 - parameters.nue)/2];
elseif strcmp(analysis.type, 'planeStrain')
    preFactor = parameters.E*(1 - parameters.nue)/(1 + parameters.nue)/(1 - 2*parameters.nue);
    C = preFactor*[1                                   parameters.nue/(1 - parameters.nue) 0
                   parameters.nue/(1 - parameters.nue) 1                                   0
                   0                                   0                                   (1 - 2*parameters.nue)/2/(1 - parameters.nue)];
else
    error('Define variable analysis.type to be either planeStress or planeStrain');
end

% Initialize output arrays
epsilon = zeros(3, numElements);
sigma = zeros(3, numElements);

%% 1. Loop over all the elements in the mesh
for iEl = 1:length(strMsh.elements(:, 1))
    %% 1i. Get the element in the mesh
    element = strMsh.elements(iEl, 2:end);
    [~, ~, element] = ...
        intersect(element', strMsh.nodes(:, 1), 'rows', 'stable');
    
    %% 1ii. Cast the element with respect to the number of its nodes
    index = isnan(element);
    element(index) = [];
    numNodesEl = length(element);
    if numNodesEl == 3
        elementType = 'linearTriangle';
    elseif numNodesEl == 4
        elementType = 'bilinearQuadrilateral';
    else
        error('Postprocessing computations for a %d-noded element has not been implemented', ...
            length(element))
    end
    
    %% 1iii. Get the nodes in the element
    nodesMsh = strMsh.nodes(element, 2:end);

    %% 1iv. Compute the basis functions and their derivatives at the Gauss Point
    if strcmp(elementType, 'linearTriangle')
        [dN, ~] = computeCST2DBasisFunctionsAndFirstDerivatives ...
            (nodesMsh(1, :), nodesMsh(2, :), nodesMsh(3, :), 0, 0);
    elseif strcmp(elementType, 'bilinearQuadrilateral')
        dN = zeros(numNodesEl, 3);
        warning('The postprocessing computations for a bilinear quadrilateral are still incomplete');
    end

    %% 1v. Compute the displacement gradient components
    uGrad = zeros(2, 2);
    for iNodeEl = 1:numNodesEl
        uGrad(1, 1) = uGrad(1, 1) + dN(iNodeEl, 2)*dHat(2*element(iNodeEl, 1) - 1, 1);
        uGrad(2, 2) = uGrad(2, 2) + dN(iNodeEl, 3)*dHat(2*element(iNodeEl, 1), 1);
        uGrad(1, 2) = uGrad(1, 2) + dN(iNodeEl, 3)*dHat(2*element(iNodeEl, 1) - 1, 1);
        uGrad(2, 1) = uGrad(2, 1) + dN(iNodeEl, 2)*dHat(2*element(iNodeEl, 1), 1);
    end

    %% 1vi. Compute the strain vector in a Voigt notation at the current element
    epsilon(:, iEl) = [uGrad(1, 1) + 1/2*(uGrad(1, 1)^2 + uGrad(2, 1)^2)
                       uGrad(2, 2) + 1/2*(uGrad(1, 2)^2 + uGrad(2, 2)^2)
                       uGrad(1, 2) + uGrad(2, 1) + uGrad(1, 1)*uGrad(1, 2) + uGrad(2, 1)*uGrad(2, 2)];

    %% 1vii. Compute the stress vector in a Voigt notation at the current element
    sigma(:, iEl) = C*epsilon(:, iEl);

    %% 1viii. Correct the shear component of the strain epsilon = [epsilonXX epsilonYY 2*epsilonXY]
    epsilon(3, iEl) = 1/2*epsilon(3, iEl);
end

end