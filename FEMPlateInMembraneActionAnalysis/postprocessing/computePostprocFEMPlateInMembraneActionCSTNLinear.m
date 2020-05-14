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
%   1iv. Create an Element Freedome Table (EFT)
%
%    1v. Compute the basis functions and their derivatives at the Gauss Point
%
%   1vi. Form the B-Operator matrix for the plate in membrane action problem
%
%  1vii. Compute the strain vector in a Voigt notation at the current element
%
% 1viii. Compute the stress vector in a Voigt notation at the current element
%
%   1ix. Correct the shear component of the strain epsilon = [epsilonXX epsilonYY 2*epsilonXY]
% <-
%
%% Function main body

%% 0. Read input

% Number of nodes in the mesh
noElements = length(strMsh.elements);

% Compute the material matrix for the given problem (The shear entry is 
% multiplied by two so that it returns the true strain field and not the 
% one needed for the computation of the internal virtual work)
if strcmp(analysis.type,'planeStress')
    preFactor = parameters.E/(1-parameters.nue^2);
    C = preFactor*[1             parameters.nue  0
                   parameters.nue 1              0
                   0              0             (1-parameters.nue)/2];
elseif strcmp(analysis.type,'planeStrain')
    preFactor = parameters.E*(1-parameters.nue)/(1+parameters.nue)/(1-2*parameters.nue);
    C = preFactor*[1                                 parameters.nue/(1-parameters.nue) 0
                   parameters.nue/(1-parameters.nue) 1                                 0
                   0                                 0                                 (1-2*parameters.nue)/2/(1-parameters.nue)];
else
    error('Define variable analysis.type to be either planeStress or planeStrain');
end

% Initialize output arrays
epsilon = zeros(3,noElements);
sigma = zeros(3,noElements);

%% 1. Loop over all the elements in the mesh
for counterEl = 1:length(strMsh.elements(:,1))
    %% 1i. Get the element in the mesh
    element = strMsh.elements(counterEl,2:end);
    
    %% 1ii. Cast the element with respect to the number of its nodes
    index = isnan(element);
    element(index) = [];
    noNodesEl = length(element);
    noDOFsEl = 2*noNodesEl;
    if noNodesEl == 3
        elementType = 'linearTriangle';
    elseif noNodesEl == 4
        elementType = 'bilinearQuadrilateral';
    else
        error('Postprocessing computations for a %d-noded element has not been implemented',length(element))
    end
    
    %% 1iii. Get the nodes in the element
    nodes = strMsh.nodes(element,2:end);
    
    %% 1iv. Create an Element Freedome Table (EFT)
    EFT = zeros(noDOFsEl,1);
    for counterEFT = 1:noNodesEl
        EFT(2*counterEFT-1) = 2*element(1,counterEFT)-1;
        EFT(2*counterEFT) = 2*element(1,counterEFT);
    end
    
    %% 1v. Compute the basis functions and their derivatives at the Gauss Point
    if strcmp(elementType,'linearTriangle')
        [dN,~] = computeCST2DBasisFunctionsAndFirstDerivatives(nodes(1,:),nodes(2,:),nodes(3,:),0,0);
    elseif strcmp(elementType,'bilinearQuadrilateral')
        dN = zeros(noNodesEl,3);
        warning('The postprocessing computations for a bilinear quadrilateral are still incomplete');
    end
    
    %% 1vi. Form the B-Operator matrix for the plate in membrane action problem
    B = zeros(3,noDOFsEl);
    for counterBOperator = 1:noNodesEl
        B(1,2*counterBOperator-1) = dN(counterBOperator,2);
        B(2,2*counterBOperator) = dN(counterBOperator,3);
        B(3,2*counterBOperator-1) = dN(counterBOperator,3);
        B(3,2*counterBOperator) = dN(counterBOperator,2);
    end
     
    %% 1vii. Compute the strain vector in a Voigt notation at the current element
    epsilon(:,counterEl) = B*dHat(EFT);
    
    %% 1viii. Compute the stress vector in a Voigt notation at the current element
    sigma(:,counterEl) = C*epsilon(:,counterEl);
    
    %% 1ix. Correct the shear component of the strain epsilon = [epsilonXX epsilonYY 2*epsilonXY]
    epsilon(3,counterEl) = epsilon(3,counterEl)/2;
end

end
