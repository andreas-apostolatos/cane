function plot_postprocResultantsIGAMembraneNLinear...
    (p,q,Xi,Eta,CP,isNURBS,parameters,DOFNumbering,xiGrid,etaGrid,dHat,graph)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Plots the reference configuration of an isogeometric Kirchhoff-Love shell
% together with the selected postprocessing resultant corresponding to
% nonlinear analysis.
%
%          Input :
%            p,q : Polynomial degrees
%         Xi,Eta : Knot vectors in xi,eta-direction
%             CP : Control point coordinates and weights of the undeformed 
%                  plate
%        isNURBS : Flag on whether the geometrical basis is NURBS or 
%                  B-Spline
%     parameters : Technical and geometrical parameters of the plate
%   DOFNumbering : The numbering of the DOFs for the given patch
%           dHat : The displacement field of the control points
% xiGrid,etaGrid : The grid points used for the plotting of the NURBS
%                  geometry
%           dHat : The displacement field of the control points
%          graph : Information on the graphics
%
%         Output :
%                  graphics
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the parametric coordinates in eta-direction
% ->
%    1i. Find the span in the eta-direction
%
%   1ii. Initialize coordinate in xi-direction
%
%  1iii. Initialize counter in xi-direction
%
%   1iv. Loop over all the parametric coordinates in xi-direction
%   ->
%        1iv.1. Find the span in xi-direction
%
%        1iv.2 Compute the IGA basis functions and possibly their derivatives
%
%        1iv.3 Compute the Cartesian image of the paratric point
%
%        1iv.4. Create the element freedom table
%
%        1iv.5. Compute the displacement field
%
%        1iv.6. Compute the covariant base vectors of the reference configuration
%
%        1iv.7. Compute the covariant base vectors of the current configuration
%
%        1iv.8. Compute the covariant metric coefficients of the current configuration
%
%        1iv.9. Compute the contravariant basis
%
%        1iv.10. Compute the local Cartesian basis
%
%        1iv.11. Compute the transformation matrix from the contravariant basis to the local Cartesian one
%
%        1iv.12. Compute the Green-Lagrange strains in the contravariant basis
%
%        1iv.13. Transform the Green-Lagrange strains in the local Cartesian basis
%
%        1iv.14. Compute the prestress values on the local Cartesian coordinate systems
%
%        1iv.15. Compute the 2nd Piola-kirchhoff stresses in the local Cartesian system
%
%        1iv.16. Compute the resultant at the evaluation point
%
%        1iv.17. Update the counter in xi-direction
%
%        1iv.18.  Update the parametric coordinate in xi-direction
%   <-
%
%    1v. Update the counter in eta-direction
%
%   1vi. Update the parametric coordinate in eta-direction
% <-
%
% 2. Plot the resultant and the knots over the domain
%
%% Function main body

%% 0. Read input

% Define a tolerance for excluding the singularity points
tolSingularity = 1e-2;

% Number of knots in u,v-direction
mxi = length(Xi);
meta = length(Eta);

% Assign a tolerance values
tol = 1e-9;

% Compute incremental steps for the resultant computation over the domain
% incremental step for eta:
deta = (Eta(meta)-Eta(1))/(xiGrid-1);

% incremental step for xi:
dxi = (Xi(mxi)-Xi(1))/(etaGrid-1);

% Initialize array of the resultant to be visualized
resultantComponent = zeros(xiGrid,etaGrid);

% Cartesian image of the parameter space
P = zeros(xiGrid,etaGrid,3);

% Number of Control Points in xi-,eta-direction
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));

% Compute the membrane material matrix (/2 at the shear component is not 
% evident for the nonlinear analysis)
materialMtxVoigt = parameters.E*parameters.t/(1-parameters.nue^2)*...
    [1              parameters.nue 0
     parameters.nue 1              0
     0              0              1-parameters.nue];
 
% Thickness of the membrane
thickness = parameters.t;
 
% Compute number of Control Points
noNodesLoc = (p + 1)*(q + 1);

% Compute number of DOFs
noDOFsLoc = 3*noNodesLoc;

% Initialize auxiliary arrays
BDisplacementsGC = zeros(3,noDOFsLoc);
EFT = zeros(1,noDOFsLoc);

% Prestress components in the local Cartesian basis
prestress = parameters.prestress;

% Compute the control point coordinates for the deformed configuration
CPd = computeDisplacedControlPointsForIGAKirchhoffLoveShell(CP,dHat);

% Initialize flags
isUndeformed = 0;

% Initialize counting index in eta-direction
counterEta = 1;  

% Initialize coordinate in eta-direction
eta = Eta(1);

%% 1. Loop over all the parametric coordinates in eta-direction
while eta <= Eta(meta) + tol
    %% 1i. Find the span in the eta-direction
    etaSpan = findKnotSpan(eta,Eta,neta);
    
    %% 1ii. Initialize coordinate in xi-direction
    xi = Xi(1);
    
    %% 1iii. Initialize counter in xi-direction
    counterXi = 1;
    
    %% 1iv. Loop over all the parametric coordinates in xi-direction
    while xi <= Xi(mxi)+tol
        %% 1iv.1. Find the span in xi-direction
        xiSpan = findKnotSpan(xi,Xi,nxi);
        
        %% 1iv.2 Compute the IGA basis functions and possibly their derivatives
        if strcmp(graph.resultant,'displacement')
            noDrvBasis = 0;
        elseif strcmp(graph.resultant,'strain') || strcmp(graph.resultant,'force')
            noDrvBasis = 1;
        end
        dR = computeIGABasisFunctionsAndDerivativesForSurface...
            (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,noDrvBasis);
        
        %% 1iv.3 Compute the Cartesian image of the paratric point
        P(counterXi,counterEta,1:3) = computeCartesianCoordinatesOfAPointOnBSplineSurface...
            (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,dR(:,1));
        
        %% 1iv.4. Create the element freedom table
        r = 1;
        for cpj = etaSpan - q:etaSpan
            for cpi = xiSpan - p:xiSpan
                EFT(r)   = DOFNumbering(cpi,cpj,1);
                EFT(r + 1) = DOFNumbering(cpi,cpj,2);
                EFT(r + 2) = DOFNumbering(cpi,cpj,3);
                r = r + 3;
            end
        end
        
        %% 1iv.5. Compute the displacement field
        if strcmp(graph.resultant,'displacement')
            k = 0;
            for c = 0:q
                for b = 0:p
                    % Update counter
                    k = k + 1;

                    % Matrix containing the basis functions
                    BDisplacementsGC(1,3*k - 2) = dR(k,1);
                    BDisplacementsGC(2,3*k - 1) = dR(k,1);
                    BDisplacementsGC(3,3*k) = dR(k,1);
                end
            end
            dispVct = BDisplacementsGC*dHat(EFT);
        end
        
        %% 1iv.6. Compute the covariant base vectors of the reference configuration
        if strcmp(graph.resultant,'strain') || strcmp(graph.resultant,'force')
            [A1,A2] = computeBaseVectorsAndDerivativesForBSplineSurface...
                (xiSpan,p,etaSpan,q,CP,0,dR);
        end
        
        %% 1iv.7. Compute the covariant base vectors of the current configuration
        if strcmp(graph.resultant,'strain') || strcmp(graph.resultant,'force')
            [a1,a2] = computeBaseVectorsAndDerivativesForBSplineSurface...
                (xiSpan,p,etaSpan,q,CPd,0,dR);
        end
        
        %% 1iv.8. Compute the covariant metric coefficients of the current configuration
        if strcmp(graph.resultant,'strain') || strcmp(graph.resultant,'force')
            aabCovariant = [a1 a2]'*[a1 a2];
        end
        
        %% 1iv.9. Compute the contravariant basis
        if strcmp(graph.resultant,'strain') || strcmp(graph.resultant,'force')
            AabCovariant = [A1 A2]'*[A1 A2];
            if det(AabCovariant) < tolSingularity
                isSingularPoint = true;
            else
                isSingularPoint = false;
            end
            AContravariant = AabCovariant\[A1 A2]';
            AContravariant = AContravariant';
        end
        
        %% 1iv.10. Compute the local Cartesian basis
        if strcmp(graph.resultant,'strain') || strcmp(graph.resultant,'force')
            eLC = computeLocalCartesianBasis4BSplineSurface...
                ([A1 A2],AContravariant);
        end
        
        %% 1iv.11. Compute the transformation matrix from the contravariant basis to the local Cartesian one
        if strcmp(graph.resultant,'strain') || strcmp(graph.resultant,'force')
            TFromContraToLC4VoigtStrain = ...
                computeTFromContra2LocalCartesian4VoigtStrainIGAKLShell...
                (eLC,AContravariant);
        end
        
        %% 1iv.12. Compute the Green-Lagrange strains in the contravariant basis
        if strcmp(graph.resultant,'strain') || strcmp(graph.resultant,'force')
            EpsilonContra = .5*[aabCovariant(1,1) - AabCovariant(1,1)
                                aabCovariant(2,2) - AabCovariant(2,2)
                                aabCovariant(1,2) - AabCovariant(1,2)];
        end
        
        %% 1iv.13. Transform the Green-Lagrange strains in the local Cartesian basis
        if strcmp(graph.resultant,'strain') || strcmp(graph.resultant,'force')
            EpsilonLC = TFromContraToLC4VoigtStrain*EpsilonContra;
        end
        
        %% 1iv.14. Compute the prestress values on the local Cartesian coordinate systems
        if strcmp(graph.resultant,'force')
            % Check if a user defined coordinate system for the prestresses 
            % is chosen
            isPrestressOverDefinedSystem = false;
            if isfield(parameters.prestress,'computeBaseVectors')
                if ~isfield(parameters.prestress,'computeParametricCoordinates')
                    error('Function handle parameters.prestress.computeParametricCoordinates has to be defined when defining the prestress over a user-defined coordinate system');
                end
                isPrestressOverDefinedSystem = true;
            end

            % Compute the convective coordinates of the surface
            if isPrestressOverDefinedSystem || isa(parameters.prestress.voigtVector,'function_handle')
                X = computeCartesianCoordinatesOfAPointOnBSplineSurface...
                    (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,dR(:,1));
                theta = prestress.computeParametricCoordinates(X);
            end

            % Compute the transformation matrix from the user defined 
            % coordinate system to the local Cartesian coordinate system if 
            % a user defined coordinate system is chosen
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

            % Transform the vector to the local Cartesian space if defined 
            % over a user defined coordinate system
%             TAndreas = [0 1 0
%                         1 0 0
%                         0 0 1];
%             pTilde = TAndreas*T2LC*pTilde;
            pTilde = T2LC*pTilde;
        end
        
        %% 1iv.15. Compute the 2nd Piola-kirchhoff stresses in the local Cartesian system
        if strcmp(graph.resultant,'force')
            NLC = thickness*pTilde +  materialMtxVoigt*EpsilonLC;
        end
        
        %% 1iv.16. Compute the resultant at the evaluation point
        if strcmp(graph.resultant,'displacement')
            if strcmp(graph.component,'x')
                resultantComponent(counterXi,counterEta) = dispVct(1,1);
            elseif strcmp(graph.component,'y')
                resultantComponent(counterXi,counterEta) = dispVct(2,1);
            elseif strcmp(graph.component,'z')
                resultantComponent(counterXi,counterEta) = dispVct(3,1);
            elseif strcmp(graph.component,'2norm')
                resultantComponent(counterXi,counterEta) = norm(dispVct);
            end
        elseif strcmp(graph.resultant,'strain')
            if isSingularPoint
                EpsilonLC = zeros(3,1);
            end
            if strcmp(graph.component,'1')
                resultantComponent(counterXi,counterEta) = EpsilonLC(1,1);
            elseif strcmp(graph.component,'2')
                resultantComponent(counterXi,counterEta) = EpsilonLC(2,1);
            elseif strcmp(graph.component,'12')
                resultantComponent(counterXi,counterEta) = EpsilonLC(3,1);
            elseif strcmp(graph.component,'1Principal')
                resultantComponent(counterXi,counterEta) = 0.5*(EpsilonLC(1,1) + ...
                    EpsilonLC(2,1) + sqrt((EpsilonLC(1,1) - EpsilonLC(2,1))^2 + ...
                    4*EpsilonLC(3,1)^2));
            elseif strcmp(graph.component,'2Principal')
                resultantComponent(counterXi,counterEta) = 0.5*(EpsilonLC(1,1) + ...
                    EpsilonLC(2,1) - sqrt((EpsilonLC(1,1) - EpsilonLC(2,1))^2 + ...
                    4*EpsilonLC(3,1)^2));
            end
        elseif strcmp(graph.resultant,'force')
            if isSingularPoint
                NLC = zeros(3,1);
            end
            if strcmp(graph.component,'1')
                resultantComponent(counterXi,counterEta) = NLC(1,1);
            elseif strcmp(graph.component,'2')
                resultantComponent(counterXi,counterEta) = NLC(2,1);
            elseif strcmp(graph.component,'12')
                resultantComponent(counterXi,counterEta) = NLC(3,1);
            elseif strcmp(graph.component,'1Principal')
                resultantComponent(counterXi,counterEta) = 0.5*(NLC(1,1)+NLC(2,1) + ...
                    sqrt((NLC(1,1)-NLC(2,1))^2 + 4*NLC(3,1)^2));
            elseif strcmp(graph.component,'2Principal')
                resultantComponent(counterXi,counterEta) = 0.5*(NLC(1,1)+NLC(2,1) - ...
                    sqrt((NLC(1,1)-NLC(2,1))^2 + 4*NLC(3,1)^2));
            end
        end

        %% 1iv.17. Update the counter in xi-direction
        counterXi = counterXi + 1;
        
        %% 1iv.18.  Update the parametric coordinate in xi-direction
        xi = xi + dxi;
    end
    
    %% 1v. Update the counter in eta-direction
    counterEta = counterEta + 1;
    
    %% 1vi. Update the parametric coordinate in eta-direction
    eta = eta + deta;
end

%% 2. Plot the resultant and the knots over the domain

% Plot the resultant over the reference configuration
% resultantComponent = zeros(xiGrid,etaGrid);
surf(P(:,:,1),P(:,:,2),P(:,:,3),resultantComponent(:,:));
hold on;

% Plot the element boundaries on the undeformed geometry
plot_knotsForBSplineSurfaceOnCartesianSpace(p,q,Xi,Eta,CP,isNURBS,isUndeformed,xiGrid,etaGrid);
hold off;

end
