function [index,minComp,maxComp] = plot_currentConfigurationAndResultants...
    (mesh,homDBC,displacement,parameters,analysis,resultant,component,graph)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Plots the current configuration of a plate in membrane action given the
% displacement field and visualizes the displacement/strain or stress field
% over the initial configuration of the plate
%
%              input :
%               mesh : Elements and nodes of the mesh
%             homDBC : Vector of the Dirichlet boundary conditions with their
%                      global numbering
%       displacement : The displacement field sorted in a vector according to its
%                      global numbering
%         parameters : The material properties of the structure
%           analysis : Analysis type (plane stress or plane strain)
%          resultant : Resultant to be visualized
%          component : component of the resultant to be visualized
%              graph : Structure on the graphics (indices etc.)
%
%             output :
%              index : The index of the current graph
%            minComp : Minmum component of the chosen component of the
%                      selected resultant in an absolute value sense
%            maxComp : Maximum component of the chosen component of the
%                      selected resultant in an absolute value sense
%
% Function layout :
%
% 0. Read input
%
% 1. Compute the displaced locations for the vertices of the triangles in the mesh
%
% 2. In the 1st window plot the deformed and the undeformed geometry
%
% 3. Plot the reference configuration
%
% 4. Plot the Current configuration
%
% 5. Plot the Dirichlet boundary conditions on the mesh
%
% 6. Assign the graph properties
%
% 7. In the 1st window plot the chosen resultant over the domain
%
% 8. Plot the mesh edges
%
% 9. Loop over all elements in the mesh
% ->
%    9i. Initialize the array of the grid points and the tensorial quantity to be visualized
%
%   9ii. Initialize local coordinates
%
%  9iii. Get element from the mesh
%
%   9iv. Get coordinates of the element vertices
%
%    9v. Get the vertices of the current element
%
%   9vi. Create the element freedom table
%
%  9vii. Get the element displacement vector
%
% 9viii. Get the moving vertices
%
%   9ix. Loop over all the sampling points
%   ->
%        9ix.1. Initialize the displacement components
%
%        9ix.2. Get the point in the interior of the line defined from Pi and Pj
%
%        9ix.3. Evaluate the CST basis functions at the evaluation point
%
%        9ix.4. Evaluate the tensorial field value on P by looping over all the products of the basis functions with the nodal solution
%
%        9ix.5. Update local coordinate by the step size
%   <-
%
%    9x. Plot the tensorial quantity over the element
% <-
%
% 10. Assign the graphics options
%
% 11. Update graph index
%
%% Function main body

%% 0. Read input

% Initialize flag on plotting
isPlotting = false;
if ~isempty(resultant)
    if ~ischar(resultant)
        error('Variable "resultant" has to be a string');
    else
        if ~strcmp(resultant,'displacement') && ~strcmp(resultant,'strain') && ~strcmp(resultant,'stress')
            error('Variable "resultant" can be one of the following; displacement, strain, stress but not %s',resultant);
        end
    end
    if ~isempty(component)
        isPlotting = true;
    end
else
    error('Variable "resultant" can not be empty');
end

% Number of nodes per element
if length(mesh.elements(1,:)) == 3
    noNodesElement = 3;
    elementType = 'triangle';
elseif length(mesh.elements(1,:)) == 4
    noNodesElement = 4;
    elementType = 'quadrilateral';
else
    error('Only triangular and quadrilateral elements are considered');
end

% Number of DoFs at the element level (depends on the element type)
isAnalysis3d = false;
if ~isempty(analysis) 
    if isfield(analysis,'type')
        if strcmp(analysis.type,'planeStress') || strcmp(analysis.type,'planeStrain')
            noDoFsElement = noNodesElement*2;
        end
    else
        noDoFsElement = noNodesElement*3;
        isAnalysis3d = true;
    end
else
    noDoFsElement = noNodesElement*3;
    isAnalysis3d = true;
end

% Compute the material matrix
if ~isempty(parameters)
    if isAnalysis3d
        thickness = parameters.t;
        materialMtxVoigt = parameters.E*thickness/(1-parameters.nue^2)*...
            [1              parameters.nue 0
             parameters.nue 1              0
             0              0              1-parameters.nue];
         isPrestress = false;
         if isfield(parameters,'prestress')
            isPrestress = true;
            prestress = parameters.prestress;
         end
    end
end

% Compute the material matrix in case strains or stresses are requested
if ~isempty(analysis)
    if isfield(analysis,'type')
        if strcmp(analysis.type,'planeStress')
            preFactor = parameters.E/(1-parameters.nue^2);
            C = preFactor*[1                parameters.nue 0
                           parameters.nue   1              0
                           0                0              (1-parameters.nue)/2];
        elseif strcmp(analysis.type,'planeStrain')
            preFactor = parameters.E*(1-parameters.nue)/(1+parameters.nue)/(1-2*parameters.nue);
            C = preFactor*[1                                 parameters.nue/(1-parameters.nue) 0
                           parameters.nue/(1-parameters.nue) 1                                 0
                           0                                 0                                (1-2*parameters.nue)/2/(1-parameters.nue)];
        end
    end
end

% Initialize output variables
minComp = inf;
maxComp = -inf;

% Define colors
if isAnalysis3d
    colorDomain = [217 218 219]/255;
else
    colorDomain = 'g';
end
% colorDomain = 'none';
% colorEdge = 'black';
% colorEdge = 'red';
colorEdge = 'none';

% Grid at each triangular element
gridXi = 1; % 5, 48
gridEta = 1; % 5, 48

% step size to both directions
stepXi = 1/(gridXi + 1);
stepEta = 1/(gridEta + 1);

%% 1. Compute the displaced locations for the vertices of the triangles in the mesh

% Initialize the array of the displaced nodes
if ~isAnalysis3d
    nodesCurrent = zeros(length(mesh.nodes),2);
else
    nodesCurrent = zeros(length(mesh.nodes),3);
end

% Initialize pseudocounter
counter = 1;

for iXi = 1:length(mesh.nodes)
    % Add the x and y components of the displacement field
    if ~isAnalysis3d
        nodesCurrent(iXi,1) = mesh.nodes(iXi,1) + displacement(2*counter - 1);
        nodesCurrent(iXi,2) = mesh.nodes(iXi,2) + displacement(2*counter);
    else
        nodesCurrent(iXi,1) = mesh.nodes(iXi,1) + displacement(3*counter - 2);
        nodesCurrent(iXi,2) = mesh.nodes(iXi,2) + displacement(3*counter - 1);
        nodesCurrent(iXi,3) = mesh.nodes(iXi,3) + displacement(3*counter);
    end
    
    % Update counter
    counter = counter + 1;
end

% Initialize figure
if isPlotting
    figure(graph.index)
end

%% 2. In the 1st window plot the deformed and the undeformed geometry
% if isPlotting
%     subplot(2,1,1);
% end
 
%% 3. Plot the Current configuration
% if isPlotting
%     if strcmp(graph.visualization.geometry,'reference_and_current') || strcmp(graph.visualization.geometry,'current');
%         patch('faces',mesh.elements,'vertices',nodesCurrent,'facecolor',colorDomain,'edgecolor',colorEdge,'LineWidth',0.01);
%         hold on;
%     end
% end

%% 4. Plot the reference configuration
% if isPlotting
%     if strcmp(graph.visualization.geometry,'reference_and_current') || strcmp(graph.visualization.geometry,'reference');
%         patch('faces',mesh.elements,'vertices',mesh.nodes,'facecolor',colorDomain,'edgecolor',colorEdge);
%     end
%     axis equal;
%     axis on;
% end
 
%% 5. Plot the Dirichlet boundary conditions on the mesh
% if isPlotting
%     if ~isempty(rb)
%         [xs,ys,zs] = createSupports(nodes_displaced,rb);
%     end
%     if ~isempty(rb)
%         hold on;
%         for k =1:length(xs(:,1))
%             plot3(xs(k,:),ys(k,:),zs(k,:),'Linewidth',2,'Color','black');
%         end
%         hold off;
%     end
% end

%% 6. Assign the graph properties
% if isPlotting
%     if ~isAnalysis3d
%         view (2);
%     end
%     camlight right;
%     lighting phong;
%     axis equal;
%     axis on;
%     xlabel('x','FontSize',14);
%     ylabel('y','FontSize',14);
%     if isAnalysis3d
%         zlabel('z','FontSize',14);
%     end
% 
%     % Title
%     if ~isAnalysis3d
%         if ~isempty(analysis)
%             if strcmp(analysis.type,'PLANE_STESS')
%                 title('Deformed configuration corresponding to plain stress analysis');
%             elseif strcmp(analysis.type,'PLANE_STRAIN')
%                 title('Deformed configuration corresponding to plain strain analysis');
%             end
%         end
%     end
% end

%% 7. In the 1st window plot the chosen resultant over the domain
% if isPlotting
%     subplot(2,1,2);
% end

%% 8. Plot the mesh edges
if isPlotting
    plot_EdgesTriangularMesh2D(mesh);
    hold on;
end

%% 9. Loop over all elements in the mesh
for iElmnts = 1:size(mesh.elements(:,1))
    %% 9i. Initialize the array of the grid points and the tensorial quantity to be visualized
    P = zeros(gridXi + 2,gridEta + 2,3);
    tensorialQuantityOverDomain = zeros(gridXi + 2,gridEta + 2);

    %% 9ii. Initialize local coordinates
    xi = 0;
    eta = 0;
    
    %% 9iii. Get element from the mesh
    element = mesh.elements(iElmnts,:);
    
    %% 9iv. Get coordinates of the element vertices
    if strcmp(elementType,'triangle')
        nodeI = element(1,1);
        nodeJ = element(1,2);
        nodeK = element(1,3);
    elseif strcmp(elementType,'quadrilateral')
        nodeI = element(1,1);
        nodeJ = element(1,2);
        nodeK = element(1,3);
        nodeL = element(1,4);
    end
    
    %% 9v. Get the vertices of the current element
    if strcmp(elementType,'triangle') || strcmp(elementType,'quadrilateral')
        Pi = mesh.nodes(nodeI,:);
        Pj = mesh.nodes(nodeJ,:);
        Pk = mesh.nodes(nodeK,:);
        pi = nodesCurrent(nodeI,:);
        pj = nodesCurrent(nodeJ,:);
        pk = nodesCurrent(nodeK,:);
    end
    if strcmp(elementType,'quadrilateral') 
        Pl = mesh.nodes(nodeL,:);
        pl = nodesCurrent(nodeL,:);
    end
    
    %% 9vi. Create the element freedom table
    if ~isAnalysis3d
        EFT = zeros(noDoFsElement,1);
        for iEta = 1:noNodesElement
            EFT(2*iEta-1,1) = 2*element(iEta)-1;
            EFT(2*iEta,1) = 2*element(iEta);
        end
    else
        EFT = zeros(noDoFsElement,1);
        for iEta = 1:noNodesElement
            EFT(3*iEta - 2,1) = 3*element(iEta) - 2;
            EFT(3*iEta - 1,1) = 3*element(iEta) - 1;
            EFT(3*iEta,1) = 3*element(iEta);
        end
    end
    
    %% 9vii. Get the element displacement vector
    displacementElement = displacement(EFT);

    %% 9viii. Get the moving vertices
    if strcmp(elementType,'triangle')
        Pi_eta = Pi;
        Pj_eta = Pj;
    elseif strcmp(elementType,'quadrilateral')
        [GPXi,~] = getGaussPointsAndWeightsOverUnitDomain(gridXi + 2);
        [GPEta,~] = getGaussPointsAndWeightsOverUnitDomain(gridEta + 2);
    end
    
    %% 9ix. Loop over all the sampling points
    for iEta = 1:gridEta + 2
        for iXi = 1:gridXi + 2
            %% 9ix.1. Initialize the displacement components
            if ~isAnalysis3d
                ux = 0;
                uy = 0;
            else
                ux = 0;
                uy = 0;
                uz = 0;
            end
            
            %% 9ix.2. Get the point in the interior of the line defined from Pi and Pj
            if strcmp(elementType,'triangle')
                P(iXi,iEta,:) = xi*Pi_eta + (1 - xi)*Pj_eta;
            elseif strcmp(elementType,'quadrilateral')
                if strcmp(resultant,'displacement')
                    [N,isInside] = computeBilinearBasisFunctions(GPXi(iXi),GPEta(iEta));
                elseif strcmp(resultant,'strain') || strcmp(resultant,'stress')
                    [dN,isInside] = computeBilinearBasisFunctionsAndFirstDerivatives(GPXi(iXi),GPEta(iEta));
                    N = dN(:,1);
                    dN = dN(:,2:3);
                end
                if ~isInside
                    error('Something went wrong when computing the basis functions in the quadrilateral');
                end
                P(iXi,iEta,:) = N(1,1)*Pi + N(2,1)*Pj + N(3,1)*Pk + N(4,1)*Pl;
            else
                error('The base element can either be triangle or quadrilateral');
            end
        
            %% 9ix.3. Evaluate the CST basis functions at the evaluation point
            if strcmp(elementType,'triangle')
                if strcmp(resultant,'displacement')
                    N = computeCST2DBasisFunctions(Pi,Pj,Pk,P(iXi,iEta,1),P(iXi,iEta,2));
                elseif strcmp(resultant,'strain') || strcmp(resultant,'stress')
                    [dN,~] = computeCST2DBasisFunctionsAndFirstDerivatives(Pi,Pj,Pk,P(iXi,iEta,1),P(iXi,iEta,2));
                end
            end
            
            %% 9ix.4. Evaluate the tensorial field value on P by looping over all the products of the basis functions with the nodal solution
            if strcmp(resultant,'displacement')
                % Compute the displacement field
                for k = 1:noNodesElement
                    if ~isAnalysis3d
                        ux = ux + N(k) * displacementElement(2*k-1);
                        uy = uy + N(k) * displacementElement(2*k);
                    else
                        ux = ux + N(k) * displacementElement(3*k - 2);
                        uy = uy + N(k) * displacementElement(3*k - 1);
                        uz = uz + N(k) * displacementElement(3*k);
                    end
                end
                
                % Assign the component variable
                if ~isempty(component)
                    if strcmp(component,'x')
                        tensorialQuantityOverDomain(iXi,iEta) = ux;
                    elseif strcmp(component,'y')
                        tensorialQuantityOverDomain(iXi,iEta) = uy;
                    elseif strcmp(component,'z')
                        if isAnalysis3d
                            tensorialQuantityOverDomain(iXi,iEta) = uz;
                        else
                            error('Z-component of the displacement is requested but no 3d analysis is assumed');
                        end
                    elseif strcmp(component,'2norm')
                        if ~isAnalysis3d
                            tensorialQuantityOverDomain(iXi,iEta) = sqrt(ux^2 + uy^2);
                        else
                            tensorialQuantityOverDomain(iXi,iEta) = sqrt(ux^2 + uy^2 + uz^2);
                        end
                    end
                    if abs(tensorialQuantityOverDomain(iXi,iEta)) < minComp
                        minComp = abs(tensorialQuantityOverDomain(iXi,iEta));
                    end
                    if abs(tensorialQuantityOverDomain(iXi,iEta)) > maxComp
                        maxComp = abs(tensorialQuantityOverDomain(iXi,iEta));
                    end
                else
                    if isAnalysis3d
                        if min([norm(ux) norm(uy) norm(uz)]) < minComp
                            minComp = min([norm(ux) norm(uy) norm(uz)]);
                        end
                        if max([norm(ux) norm(uy) norm(uz)]) > maxComp
                            maxComp = max([norm(ux) norm(uy) norm(uz)]);
                        end
                    else
                        if min([norm(ux) norm(uy)]) < minComp
                            minComp = min([norm(ux) norm(uy)]);
                        end
                        if max([norm(ux) norm(uy)]) > maxComp
                            maxComp = max([norm(ux) norm(uy)]);
                        end
                    end
                end
            elseif strcmp(resultant,'strain')
                if ~isAnalysis3d
                    % Compute the B-Operator Matrix
                    B = [dN(1,2) 0       dN(2,2) 0       dN(3,2) 0
                         0       dN(1,3) 0       dN(2,3) 0       dN(3,3)
                         dN(1,3) dN(1,2) dN(2,3) dN(2,2) dN(3,3) dN(3,2)];

                    % Compute the strain field in Voigt notation
                    strain = B*displacementElement;
                    if strcmp(component,'x')
                        tensorialQuantityOverDomain(iXi,iEta) = strain(1);
                    elseif strcmp(component,'y')
                        tensorialQuantityOverDomain(iXi,iEta) = strain(2);
                    elseif strcmp(component,'xy')
                        tensorialQuantityOverDomain(iXi,iEta) = strain(3);
                    elseif strcmp(component,'1Principal')
                        tensorialQuantityOverDomain(iXi,iEta) = 0.5*(strain(1,1) + ...
                            strain(2,1) - sqrt((strain(1,1) - strain(2,1))^2 + ...
                            strain(3,1)^2));
                    elseif strcmp(component,'2Principal')
                        tensorialQuantityOverDomain(iXi,iEta) = 0.5*(strain(1,1) + ...
                            strain(2,1) - sqrt((strain(1,1) - strain(2,1))^2 + ...
                            strain(3,1)^2));
                    end
                else
                    % Compute the base vectors of the reference
                    % configuration
                    A1 = Pj' - Pi';
                    A2 = Pk' - Pi';
                    
                    % Compute the base vectors of the current
                    % configuration
                    a1 = pj' - pi';
                    a2 = pk' - pi';
                    
                    % Compute the covariant metric coefficients of the
                    % reference configuration
                    AabCovariant = [A1 A2]'*[A1 A2];
                    
                    % Compute the covariant metric coefficients of the
                    % current configuration
                    aabCovariant = [a1 a2]'*[a1 a2];
                    
                    % Compute the contravariant basis of the reference 
                    % configuration
                    AContravariant = AabCovariant\[A1 A2]';
                    AContravariant = AContravariant';
                    
                    % Compute the local Cartesian basis
                    eLC = computeLocalCartesianBasis4BSplineSurface...
                        ([A1 A2],AContravariant);
                    
                    % Compute the transformation matrix from the 
                    % contravariant basis to the local Cartesian one
                    TFromContraToLC4VoigtStrain = zeros(3,3);
                    TFromContraToLC4VoigtStrain(1,1) = (AContravariant(:,1)'*eLC(:,1))^2;
                    TFromContraToLC4VoigtStrain(1,2) = (AContravariant(:,2)'*eLC(:,1))^2;
                    TFromContraToLC4VoigtStrain(1,3) = 2*(AContravariant(:,1)'*eLC(:,1))*(AContravariant(:,2)'*eLC(:,1));
                    TFromContraToLC4VoigtStrain(2,1) = (AContravariant(:,1)'*eLC(:,2))^2;
                    TFromContraToLC4VoigtStrain(2,2) = (AContravariant(:,2)'*eLC(:,2))^2;
                    TFromContraToLC4VoigtStrain(2,3) = 2*(AContravariant(:,1)'*eLC(:,2))*(AContravariant(:,2)'*eLC(:,2));
                    TFromContraToLC4VoigtStrain(3,1) = (AContravariant(:,1)'*eLC(:,1))*(AContravariant(:,1)'*eLC(:,2));
                    TFromContraToLC4VoigtStrain(3,2) = (AContravariant(:,2)'*eLC(:,1))*(AContravariant(:,2)'*eLC(:,2));
                    TFromContraToLC4VoigtStrain(3,3) = (AContravariant(:,1)'*eLC(:,1))*(AContravariant(:,2)'*eLC(:,2)) + ...
                        (AContravariant(:,2)'*eLC(:,1))*(AContravariant(:,1)'*eLC(:,2));
                    
                    % Compute the Green-Lagrange elastic strains in the 
                    % contravariant basis
                    EpsilonContra = .5*[aabCovariant(1,1) - AabCovariant(1,1)
                                        aabCovariant(2,2) - AabCovariant(2,2)
                                        aabCovariant(1,2) - AabCovariant(1,2)];
                                    
                    % Transform the Green-Lagrange elastic strains in the 
                    % local Cartesian basis
                    EpsilonLC = TFromContraToLC4VoigtStrain*EpsilonContra;
                    
                    % Compute the prestress values on the local Cartesian 
                    % coordinate systems
                    
                    % Check if a user defined coordinate system for the 
                    % prestresses is chosen
                    isPrestressOverDefinedSystem = false;
                    if isfield(prestress,'computeBaseVectors')
                        if ~isfield(prestress,'computeParametricCoordinates')
                            error('Function handle prestress.computeParametricCoordinates has to be defined when defining the prestress over a user-defined coordinate system');
                        end
                        isPrestressOverDefinedSystem = true;
                    end
                    
                    % Compute the convective coordinates of the surface
                    if isPrestressOverDefinedSystem || isa(prestress.voigtVector,'function_handle')
                        theta = prestress.computeParametricCoordinates(squeeze(P(iXi,iEta,:)));
                    end
                    
                    % Compute the transformation matrix from the user 
                    % defined coordinate system to the local Cartesian 
                    % coordinate system if a user defined coordinate system 
                    % is chosen
                    if isPrestressOverDefinedSystem
                        prestressBaseVct = prestress.computeBaseVectors(theta(1,1),theta(2,1));
                        T2LC = computeT2LocalCartesianBasis(prestressBaseVct,eLC);
                    else
                        T2LC = [1 0 0
                                0 1 0
                                0 0 1];
                    end
                    
                    % Compute the prestress values
                    if isa(prestress.voigtVector,'function_handle')
                        pTilde = prestress.voigtVector(theta);
                    else
                        pTilde = prestress.voigtVector;
                    end
                    
                    % Transform the vector to the local Cartesian space if 
                    % defined over a user defined coordinate system
%                     TAndreas = [0 1 0
%                                 1 0 0
%                                 0 0 1];
%                     pTilde = TAndreas*T2LC*pTilde;
                    pTilde = T2LC*pTilde;
                    
                    % Compute the total strains (elastic plus from 
                    % prestress) in the local Cartesian basis
                    EpsilonLC = EpsilonLC + materialMtxVoigt\(pTilde*thickness);
                    
                    % Save the chosen strain component at the evaluation
                    % point
                    if ~isempty(component)
                        if strcmp(component,'x')
                            tensorialQuantityOverDomain(iXi,iEta) = EpsilonLC(1,1);
                        elseif strcmp(component,'y')
                            tensorialQuantityOverDomain(iXi,iEta) = EpsilonLC(2,1);
                        elseif strcmp(component,'xy')
                            tensorialQuantityOverDomain(iXi,iEta) = EpsilonLC(3,1);
                        elseif strcmp(component,'2norm')
                            tensorialQuantityOverDomain(iXi,iEta) = sqrt(EpsilonLC(1,1)^2 + EpsilonLC(2,1)^2 + 2*EpsilonLC(3,1)^2);
                        elseif strcmp(component,'1Principal')
                            tensorialQuantityOverDomain(iXi,iEta) = 0.5*(EpsilonLC(1,1) + ...
                                EpsilonLC(2,1) + sqrt((EpsilonLC(1,1) - EpsilonLC(2,1))^2 + ...
                                4*EpsilonLC(3,1)^2));
                        elseif strcmp(component,'2Principal')
                            tensorialQuantityOverDomain(iXi,iEta) = 0.5*(EpsilonLC(1,1) + ...
                                EpsilonLC(2,1) - sqrt((EpsilonLC(1,1) - EpsilonLC(2,1))^2 + ...
                                4*EpsilonLC(3,1)^2));
                        else
                            error('Choose correct strain component to visualize');
                        end
                        if abs(tensorialQuantityOverDomain(iXi,iEta)) < minComp
                            minComp = abs(tensorialQuantityOverDomain(iXi,iEta));
                        end
                        if abs(tensorialQuantityOverDomain(iXi,iEta)) > maxComp
                            maxComp = abs(tensorialQuantityOverDomain(iXi,iEta));
                        end
                    else
                        if min(abs(EpsilonLC)) < minComp
                            minComp = min(abs(EpsilonLC));
                        end
                        if max(abs(EpsilonLC)) > maxComp
                            maxComp = max(abs(EpsilonLC));
                        end
                    end
                end
            elseif strcmp(resultant,'stress')
                if ~isAnalysis3d
                    % Compute the B-Operator Matrix    
                    B = [dN(1,2) 0       dN(2,2) 0       dN(3,2) 0
                         0       dN(1,3) 0       dN(2,3) 0       dN(3,3)
                         dN(1,3) dN(1,2) dN(2,3) dN(2,2) dN(3,3) dN(3,2)];

                    % Compute the strain field in Voigt notation
                    strain = B*displacementElement;

                    % compute the stress field in Voigt notation
                    stress = C*strain;
                    
                    % Save the chosen force component at the evaluation
                    % point
                    if strcmp(component,'x')
                        tensorialQuantityOverDomain(iXi,iEta) = stress(1);
                    elseif strcmp(component,'y')
                        tensorialQuantityOverDomain(iXi,iEta) = stress(2);
                    elseif strcmp(component,'xy')
                        tensorialQuantityOverDomain(iXi,iEta) = stress(3);
                    elseif strcmp(component,'1Principal')
                        tensorialQuantityOverDomain(iXi,iEta) = 0.5*(stress(1,1) + ...
                            stress(2,1) + sqrt((stress(1,1)-stress(2,1))^2 + 4*stress(3,1)^2));
                    elseif strcmp(component,'2Principal')
                        tensorialQuantityOverDomain(iXi,iEta) = 0.5*(stress(1,1) + ...
                            stress(2,1) - sqrt((stress(1,1)-stress(2,1))^2 + 4*stress(3,1)^2));
                    else
                        error('Choose correct force component to visualize');
                    end
                else
                    % Compute the base vectors of the reference
                    % and the current configuration
                    if strcmp(elementType,'triangle')
                        A1 = Pj' - Pi';
                        A2 = Pk' - Pi';
                        a1 = pj' - pi';
                        a2 = pk' - pi';
                    elseif strcmp(elementType,'quadrilateral')
                        A1 = dN(1,1)*Pi' + dN(2,1)*Pj' + dN(3,1)*Pk' + dN(4,1)*Pl';
                        A2 = dN(1,2)*Pi' + dN(2,2)*Pj' + dN(3,2)*Pk' + dN(4,2)*Pl';
                        a1 = dN(1,1)*pi' + dN(2,1)*pj' + dN(3,1)*pk' + dN(4,1)*pl';
                        a2 = dN(1,2)*pi' + dN(2,2)*pj' + dN(3,2)*pk' + dN(4,2)*pl';
                    end
                    
                    % Compute the covariant metric coefficients of the
                    % reference configuration
                    AabCovariant = [A1 A2]'*[A1 A2];
                    
                    % Compute the covariant metric coefficients of the
                    % current configuration
                    aabCovariant = [a1 a2]'*[a1 a2];
                    
                    % Compute the contravariant basis of the reference 
                    % configuration
                    AContravariant = AabCovariant\[A1 A2]';
                    AContravariant = AContravariant';
                    
                    % Compute the local Cartesian basis
                    eLC = computeLocalCartesianBasis4BSplineSurface...
                        ([A1 A2],AContravariant);
                    
                    % Compute the transformation matrix from the 
                    % contravariant basis to the local Cartesian one
                    TFromContraToLC4VoigtStrain = zeros(3,3);
                    TFromContraToLC4VoigtStrain(1,1) = (AContravariant(:,1)'*eLC(:,1))^2;
                    TFromContraToLC4VoigtStrain(1,2) = (AContravariant(:,2)'*eLC(:,1))^2;
                    TFromContraToLC4VoigtStrain(1,3) = 2*(AContravariant(:,1)'*eLC(:,1))*(AContravariant(:,2)'*eLC(:,1));
                    TFromContraToLC4VoigtStrain(2,1) = (AContravariant(:,1)'*eLC(:,2))^2;
                    TFromContraToLC4VoigtStrain(2,2) = (AContravariant(:,2)'*eLC(:,2))^2;
                    TFromContraToLC4VoigtStrain(2,3) = 2*(AContravariant(:,1)'*eLC(:,2))*(AContravariant(:,2)'*eLC(:,2));
                    TFromContraToLC4VoigtStrain(3,1) = (AContravariant(:,1)'*eLC(:,1))*(AContravariant(:,1)'*eLC(:,2));
                    TFromContraToLC4VoigtStrain(3,2) = (AContravariant(:,2)'*eLC(:,1))*(AContravariant(:,2)'*eLC(:,2));
                    TFromContraToLC4VoigtStrain(3,3) = (AContravariant(:,1)'*eLC(:,1))*(AContravariant(:,2)'*eLC(:,2)) + ...
                        (AContravariant(:,2)'*eLC(:,1))*(AContravariant(:,1)'*eLC(:,2));
                    
                    % Compute the Green-Lagrange strains in the 
                    % contravariant basis
                    EpsilonContra = .5*[aabCovariant(1,1) - AabCovariant(1,1)
                                        aabCovariant(2,2) - AabCovariant(2,2)
                                        aabCovariant(1,2) - AabCovariant(1,2)];
                                    
                    % Transform the Green-Lagrange strains in the local 
                    % Cartesian basis
                    EpsilonLC = TFromContraToLC4VoigtStrain*EpsilonContra;
                    
                    % Compute the prestress values on the local Cartesian 
                    % coordinate systems
                    
                    % Check if a user defined coordinate system for the 
                    % prestresses is chosen
                    isPrestressOverDefinedSystem = false;
                    if isfield(prestress,'computeBaseVectors')
                        if ~isfield(prestress,'computeParametricCoordinates')
                            error('Function handle prestress.computeParametricCoordinates has to be defined when defining the prestress over a user-defined coordinate system');
                        end
                        isPrestressOverDefinedSystem = true;
                    end
                    
                    % Compute the convective coordinates of the surface
                    if isPrestressOverDefinedSystem || isa(prestress.voigtVector,'function_handle')
                        theta = prestress.computeParametricCoordinates(squeeze(P(iXi,iEta,:)));
                    end
                    
                    % Compute the transformation matrix from the user 
                    % defined coordinate system to the local Cartesian 
                    % coordinate system if a user defined coordinate system 
                    % is chosen
                    if isPrestressOverDefinedSystem
                        prestressBaseVct = prestress.computeBaseVectors(theta(1,1),theta(2,1));
                        T2LC = computeT2LocalCartesianBasis(prestressBaseVct,eLC);
                    else
                        T2LC = [1 0 0
                                0 1 0
                                0 0 1];
                    end
                    
                    % Compute the prestress values
                    if isa(prestress.voigtVector,'function_handle')
                        pTilde = prestress.voigtVector(theta);
                    else
                        pTilde = prestress.voigtVector;
                    end
                    
                    % Transform the vector to the local Cartesian space if 
                    % defined over a user defined coordinate system
%                     TAndreas = [0 1 0
%                                 1 0 0
%                                 0 0 1];
%                     pTilde = TAndreas*T2LC*pTilde;
                    pTilde = T2LC*pTilde;
                    
                    % Compute the force field resulting from the elastic
                    % strains
                    NLC = thickness*pTilde + materialMtxVoigt*EpsilonLC;
                    
                    % Save the chosen force component at the evaluation
                    % point
                    if ~isempty(component)
                        if strcmp(component,'x')
                            tensorialQuantityOverDomain(iXi,iEta) = NLC(1,1);
                        elseif strcmp(component,'y')
                            tensorialQuantityOverDomain(iXi,iEta) = NLC(2,1);
                        elseif strcmp(component,'xy')
                            tensorialQuantityOverDomain(iXi,iEta) = NLC(3,1);
                        elseif strcmp(component,'2norm')
                            tensorialQuantityOverDomain(iXi,iEta) = sqrt(NLC(1,1)^2 + NLC(2,1)^2 + 2*NLC(3,1)^2);
                        elseif strcmp(component,'1Principal')
                            tensorialQuantityOverDomain(iXi,iEta) = 0.5*(NLC(1,1) + ...
                                NLC(2,1) + sqrt((NLC(1,1) - NLC(2,1))^2 + 4*NLC(3,1)^2));
                        elseif strcmp(component,'2Principal')
                            tensorialQuantityOverDomain(iXi,iEta) = 0.5*(NLC(1,1) + ...
                                NLC(2,1) - sqrt((NLC(1,1) - NLC(2,1))^2 + 4*NLC(3,1)^2));
                        else
                            error('Choose correct force component to visualize');
                        end
                        if abs(tensorialQuantityOverDomain(iXi,iEta)) < minComp
                            minComp = abs(tensorialQuantityOverDomain(iXi,iEta));
                        end
                        if abs(tensorialQuantityOverDomain(iXi,iEta)) > maxComp
                            maxComp = abs(tensorialQuantityOverDomain(iXi,iEta));
                        end
                    else
                        if min(abs(NLC)) < minComp
                            minComp = min(abs(NLC));
                        end
                        if max(abs(NLC)) > maxComp
                            maxComp = max(abs(NLC));
                        end
                    end
                end
            end
                  
            %% 9ix.5. Update local coordinate by the step size
            xi = xi + stepXi;
        end
        xi = 0;
        eta = eta + stepEta;
        Pi_eta = eta*Pk + (1-eta)*Pi;
        Pj_eta = eta*Pk + (1-eta)*Pj;
    end

    %% 9x. Plot the tensorial quantity over the element
    if isPlotting
        surf(P(:,:,1),P(:,:,2),P(:,:,3),tensorialQuantityOverDomain(:,:,1),'EdgeColor','none');
        hold on;
    end
end

%% 10. Assign the graphics options
if isPlotting
    colormap('jet');
    colorbar;
    hold on;
    view (2);
    axis equal;
    axis on;
    grid on;
    xlabel('x','FontSize',14);
    ylabel('y','FontSize',14);
    if strcmp(resultant,'displacement')
        if strcmp(component,'x')
            title('x-component of the displacement field u_x');
        elseif strcmp(component,'y')
            title('y-component of the displacement field u_y');
        elseif strcmp(component,'z')
            title('y-component of the displacement field u_z');
        elseif strcmp(component,'2norm')
            title('2-norm of the displacement field ||u||_2');
        end
    elseif strcmp(resultant,'strain')
        if strcmp(component,'x')
            title('normal xx-component of the strain field epsilon_{xx}');
        elseif strcmp(component,'y')
            title('normal yy-component of the strain field epsilon_{yy}');
        elseif strcmp(component,'xy')
            title('shear xy-component of the strain field epsilon_{xy}');
        elseif strcmp(component,'1Principal')
            title('1st principal component of the strain field epsilon_1');
        elseif strcmp(component,'2Principal')
            title('2nd principal component of the strain field epsilon_2');
        end
    elseif strcmp(resultant,'stress')
        if strcmp(component,'x')
            title('normal xx-component of the stress field sigma_{xx}');
        elseif strcmp(component,'y')
            title('normal yy-component of the stress field sigma_{yy}');
        elseif strcmp(component,'xy')
            title('shear xy-component of the stress field sigma_{xy}');
        elseif strcmp(component,'1Principal')
            title('1st principal component of the stress field sigma_1');
        elseif strcmp(component,'2Principal')
            title('2nd principal component of the stress field sigma_2');
        end
    end
    hold off;
end

%% 11. Update graph index
if isPlotting
    index = graph.index + 1;
else
    index = graph.index;
end

end
