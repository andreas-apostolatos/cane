function plot_postprocResultantsIGAKirchhoffLoveShell...
    (p,q,Xi,Eta,CP,isNURBS,parameters,xiGrid,etaGrid,dHat,graph)
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
% together with the selected postprocessing resultant.
%
%          Input :
%            p,q : Polynomial degrees
%         Xi,Eta : Knot vectors in xi,eta-direction
%             CP : Control point coordinates and weights of the undeformed 
%                  plate
%        isNURBS : Flag on whether the geometrical basis is NURBS or 
%                  B-Spline
%     parameters : Technical and geometrical parameters of the plate
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
% 1. Compute the resultant array to be used for the visualization
%
% 2. Visualize the resultant and the knots over the domain
%
%% Function main body

%% 0. Read input

% Number of knots in u,v-direction
mxi = length(Xi);
meta = length(Eta);

% Assign a tolerance value
tol = 10e-10;

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

% Compute the membrane material matrix
Dm = parameters.E*parameters.t/(1-parameters.nue^2)*...
    [1              parameters.nue 0
	 parameters.nue 1              0
     0              0              (1-parameters.nue)/2];
                                                 
% Compute the bending material matrix
Db = parameters.E*parameters.t^3/(12*(1-parameters.nue^2))*...
    [1              parameters.nue 0
     parameters.nue 1              0 
     0              0              (1-parameters.nue)/2];
 
% Make a DOF numbering for the given patch
dHatElem = zeros(mxi-p-1,meta-q-1,3*(p+1)*(q+1));
for etaSpan = (q+1):(meta-q-1)
    for xiSpan = (p+1):(mxi-p-1)
        xiCounter = 1; 
        for c = etaSpan-q-1:etaSpan-1 
            for b = xiSpan-p:xiSpan
                dHatElem(xiSpan,etaSpan,xiCounter) = dHat(3*(c*nxi+b)-2);
                dHatElem(xiSpan,etaSpan,xiCounter + 1) = dHat(3*(c*nxi+b)-1);
                dHatElem(xiSpan,etaSpan,xiCounter + 2) = dHat(3*(c*nxi+b));
                
                % Update counter
                xiCounter = xiCounter + 3;
            end
        end
    end
end

% Compute the control point coordinates for the deformed configuration
CPd = computeDisplacedControlPointsForIGAKirchhoffLoveShell(CP,dHat);

% Initialize flags
isUndeformed = 0;
isDeformed = 0;

%% 1. Compute the resultant array to be used for the visualization

% counting index in eta-direction
etaCounter = 1;  

% Initialize coordinate in eta-direction
eta = Eta(1);

% Loop over all the parametric coordinates in eta-direction
while eta <= Eta(meta)+tol
    % Find the span in the eta-direction
    etaSpan = findKnotSpan(eta,Eta,neta);
    
    % Initialize coordinate in xi-direction
    xi = Xi(1);
    
    % Initialize counter in xi-direction
    xiCounter = 1;
    
    % Loop over all the parametric coordinates in xi-direction
    while xi <= Xi(mxi)+tol
        % Find the span in xi-direction
        xiSpan = findKnotSpan(xi,Xi,nxi);
        
        % Compute the IGA basis functions and possibly their derivatives
        if strcmp(graph.resultant,'displacement')
            nDrv = 0;
        elseif strcmp(graph.resultant,'rotation') || strcmp(graph.resultant,'strain') ...
                || strcmp(graph.resultant,'force') || strcmp(graph.resultant,'curvature') || ...
                strcmp(graph.resultant,'moment')
            nDrv = 2;
        elseif strcmp(graph.resultant,'shearForce')
            nDrv = 3;
        end
        dR = computeIGABasisFunctionsAndDerivativesForSurface(xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,nDrv);
        
        % Compute necessary metrics
        if strcmp(graph.resultant,'rotation') || strcmp(graph.resultant,'strain') ||...
                strcmp(graph.resultant,'force') || strcmp(graph.resultant,'curvature') || ...
                strcmp(graph.resultant,'moment') ||...
                strcmp(graph.resultant,'shearForce')
            % Compute the base vectors and possibly their first derivatives
            if strcmp(graph.resultant,'strain') ||...
                    strcmp(graph.resultant,'force')
                nDrvBaseVct = 0;
                indexBaseVct = [1 2 4];
            elseif strcmp(graph.resultant,'rotation') || ...
                    strcmp(graph.resultant,'curvature') || ...
                    strcmp(graph.resultant,'moment')
                nDrvBaseVct = 1;
                indexBaseVct = [1 2 3 4 5 6];
                if strcmp(graph.resultant,'rotation')
                   BDisplacementsGC = zeros(3,(p+1)*(q+1));
                   dRdxi = zeros(3,3*(p+1)*(q+1));
                   dRdeta = zeros(3,3*(p+1)*(q+1));
                end
            elseif strcmp(graph.resultant,'shearForce')
                nDrvBaseVct = 2;
                indexBaseVct = [1 2 3 4 5 6 7 8 9 10];
            end
            [dA1,dA2] = computeBaseVectorsAndDerivativesForBSplineSurface(xiSpan,p,etaSpan,q,CP,nDrvBaseVct,dR(:,indexBaseVct));
            
            % Compute the surface normal (third covariant base vector not 
            % normalized) for the reference configuration
            A3Tilde = cross(dA1(:,1),dA2(:,1));
            
            % Compute the surface normal vector
            A3 = A3Tilde/norm(A3Tilde);
        end
        
        % Compute the Cartesian image of the paratric point
        P(xiCounter,etaCounter,1:3) = computeCartesianCoordinatesOfAPointOnBSplineSurface(xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,dR(:,1));
        
        % Get the actual Control Point displacement vector for the current
        % knot span
        dHatActual(:,1) = dHatElem(xiSpan,etaSpan,:);
        
        % Compute the resultant component
        if strcmp(graph.resultant,'displacement')
            % Compute the displacement field at the given parametric
            % location
            d = computePostprocDisplacementIGAKirchhoffLoveShell(p,q,dR(:,1),dHatActual);
            
            % Decide upon the the visualization resultant
            if strcmp(graph.component,'x')
                resultantComponent(xiCounter,etaCounter) = d(1,1);
            elseif strcmp(graph.component,'y')
                resultantComponent(xiCounter,etaCounter) = d(2,1);
            elseif strcmp(graph.component,'z')
                resultantComponent(xiCounter,etaCounter) = d(3,1);
            elseif strcmp(graph.component,'2norm')
                resultantComponent(xiCounter,etaCounter) = norm(d);
            end
        elseif strcmp(graph.resultant,'strain')
            % Compute the B-operator matrix for the strain field
            BStrain = computeBOperatorMatrix4StrainIGAKirchhoffLoveShellLinear...
                (p,q,dR,[dA1(:,1) dA2(:,1)]);
            
            % Compute the strain tensor in a Voigt notation
            epsilonVoigt = BStrain*dHatActual;
            
            % Decide upon the the visualization resultant
            if strcmp(graph.component,'1')
                resultantComponent(xiCounter,etaCounter) = epsilonVoigt(1,1);
            elseif strcmp(graph.component,'2')
                resultantComponent(xiCounter,etaCounter) = epsilonVoigt(2,1);
            elseif strcmp(graph.component,'12')
                resultantComponent(xiCounter,etaCounter) = epsilonVoigt(3,1);
            elseif strcmp(graph.component,'1Principal')
                resultantComponent(xiCounter,etaCounter) = 0.5*(epsilonVoigt(1,1)+epsilonVoigt(2,1) + sqrt((epsilonVoigt(1,1)-epsilonVoigt(2,1))^2 + epsilonVoigt(3,1)^2));
            elseif strcmp(graph.component,'2Principal')
                resultantComponent(xiCounter,etaCounter) = 0.5*(epsilonVoigt(1,1)+epsilonVoigt(2,1) - sqrt((epsilonVoigt(1,1)-epsilonVoigt(2,1))^2 + epsilonVoigt(3,1)^2));
            end
        elseif strcmp(graph.resultant,'rotation')
            % Compute the derivatives of the surface normal vector
            dA3 = computeParametricDrvsSurfaceNormalOnBSplineSurface...
                ([dA1(:,1) dA2(:,1)],[dA1(:,2) dA2(:,2) dA1(:,3)],...
                A3,norm(A3Tilde));
            
            % Compute the B-operator matrix for the displacement field
            k = 1;
            for c = 0:q
                for b = 0:p
                    % Matrix containing the basis functions
                    BDisplacementsGC(1,3*k-2) = dR(k,1);
                    BDisplacementsGC(2,3*k-1) = dR(k,1);
                    BDisplacementsGC(3,3*k) = dR(k,1);
                    
                    % Matrix containing the derivatives of the basis functions
                    % With respect to xi:
                    dRdxi(1,3*k-2) = dR(k,2);
                    dRdxi(2,3*k-1) = dR(k,2);
                    dRdxi(3,3*k) = dR(k,2);
                    
                    % With respect to eta:
                    dRdeta(1,3*k-2) = dR(k,4);
                    dRdeta(2,3*k-1) = dR(k,4);
                    dRdeta(3,3*k) = dR(k,4);
                    
                    % Update counter
                    k = k + 1;
                end
            end
            
            % Compute the contravariant basis
            AabCov = [dA1(:,1) dA2(:,1)]'*[dA1(:,1) dA2(:,1)];
            AContravariant = (AabCov\[dA1(:,1) dA2(:,1)]')';
            
            % Compute the curvature coefficients
            BV = [dA1(:,2) dA2(:,2) dA1(:,3)]'*A3;
            
            % Compute the B-operator matrices for the rotation vector
            [BOmega1,BOmega2,~,~] = computePostprocBOperatorMatrix4RotationsIGAKirchhoffLoveShell...
                (BDisplacementsGC,dRdxi,dRdeta,A3Tilde,dA3,AContravariant,BV);
            
            % Compute the rotation vector in the covariant basis
            rotationVctCov = zeros(2,1);
            rotationVctCov(1,1) = BOmega1*dHatActual;
            rotationVctCov(2,1) = BOmega2*dHatActual;

            % Decide upon the the visualization resultant
            if strcmp(graph.component,'1')
                resultantComponent(xiCounter,etaCounter) = rotationVctCov(1,1);
            elseif strcmp(graph.component,'2')
                resultantComponent(xiCounter,etaCounter) = rotationVctCov(2,1);
            elseif strcmp(graph.component,'2norm')
                resultantComponent(xiCounter,etaCounter) = norm(rotationVctCov);
            end
        elseif strcmp(graph.resultant,'force')
            % Compute the B-operator matrix for the strain field
            BStrain = computeBOperatorMatrix4StrainIGAKirchhoffLoveShellLinear...
                (p,q,dR,[dA1(:,1) dA2(:,1)]);
            
            % Compute the strain tensor in a Voigt notation
            epsilonVoigt = BStrain*dHatActual;
            
            % Compute the force tensor in a Voigt notation
            forceVoigt = Dm*epsilonVoigt;
            
            % Decide upon the the visualization resultant
            if strcmp(graph.component,'1')
                resultantComponent(xiCounter,etaCounter) = forceVoigt(1,1);
            elseif strcmp(graph.component,'2')
                resultantComponent(xiCounter,etaCounter) = forceVoigt(2,1);
            elseif strcmp(graph.component,'12')
                resultantComponent(xiCounter,etaCounter) = forceVoigt(3,1);
            elseif strcmp(graph.component,'1Principal')
                resultantComponent(xiCounter,etaCounter) = 0.5*(forceVoigt(1,1)+forceVoigt(2,1) + sqrt((forceVoigt(1,1)-forceVoigt(2,1))^2 + 4*forceVoigt(3,1)^2));
            elseif strcmp(graph.component,'2Principal')
                resultantComponent(xiCounter,etaCounter) = 0.5*(forceVoigt(1,1)+forceVoigt(2,1) - sqrt((forceVoigt(1,1)-forceVoigt(2,1))^2 + 4*forceVoigt(3,1)^2));
            end
        elseif strcmp(graph.resultant,'curvature')
            % Compute the B-operator matrix for the strain field
            BCurvature = computeBOperatorMatrix4CurvatureIGAKirchhoffLoveShellLinear(p,q,dR,[dA1(:,1) dA2(:,1)],[dA1(:,2) dA2(:,2) dA1(:,3)],A3Tilde);
            
            % Compute the strain tensor in a Voigt notation
            curvatureVoigt = BCurvature*dHatActual;
            
            % Decide upon the the visualization resultant
            if strcmp(graph.component,'1')
                resultantComponent(xiCounter,etaCounter) = curvatureVoigt(1,1);
            elseif strcmp(graph.component,'2')
                resultantComponent(xiCounter,etaCounter) = curvatureVoigt(2,1);
            elseif strcmp(graph.component,'12')
                resultantComponent(xiCounter,etaCounter) = curvatureVoigt(3,1);
            elseif strcmp(graph.component,'1Principal')
                resultantComponent(xiCounter,etaCounter) = 0.5*(curvatureVoigt(1,1)+curvatureVoigt(2,1) + sqrt((curvatureVoigt(1,1)-curvatureVoigt(2,1))^2 + curvatureVoigt(3,1)^2));
            elseif strcmp(graph.component,'2Principal')
                resultantComponent(xiCounter,etaCounter) = 0.5*(curvatureVoigt(1,1)+curvatureVoigt(2,1) - sqrt((curvatureVoigt(1,1)-curvatureVoigt(2,1))^2 + curvatureVoigt(3,1)^2));
            end
        elseif strcmp(graph.resultant,'moment')
            % Compute the B-operator matrix for the change in the curvature 
            % field
            BCurvature = computeBOperatorMatrix4CurvatureIGAKirchhoffLoveShellLinear(p,q,dR,[dA1(:,1) dA2(:,1)],[dA1(:,2) dA2(:,2) dA1(:,3)],A3Tilde);
            
            % Compute the strain tensor in a Voigt notation
            curvatureVoigt = BCurvature*dHatActual;
            
            % Compute the force tensor in a Voigt notation
            momentVoigt = Db*curvatureVoigt;
            
            % Decide upon the the visualization resultant
            if strcmp(graph.component,'1')
                resultantComponent(xiCounter,etaCounter) = momentVoigt(1,1);
            elseif strcmp(graph.component,'2')
                resultantComponent(xiCounter,etaCounter) = momentVoigt(2,1);
            elseif strcmp(graph.component,'12')
                resultantComponent(xiCounter,etaCounter) = momentVoigt(3,1);
            elseif strcmp(graph.component,'1Principal')
                resultantComponent(xiCounter,etaCounter) = 0.5*(momentVoigt(1,1)+momentVoigt(2,1) + sqrt((momentVoigt(1,1)-momentVoigt(2,1))^2 + 4*momentVoigt(3,1)^2));
            elseif strcmp(graph.component,'2Principal')
                resultantComponent(xiCounter,etaCounter) = 0.5*(momentVoigt(1,1)+momentVoigt(2,1) - sqrt((momentVoigt(1,1)-momentVoigt(2,1))^2 + 4*momentVoigt(3,1)^2));
            end
        elseif strcmp(graph.resultant,'shearForce')
            % Compute the B-operator matrix for the shear forces in the
            % local Cartesian basis
            BShearForces = computeBOperatorMatrix4ShearForcesIGAKLShellAndreasLinear...
                (p,q,dR,dA1,dA2,A3Tilde,Db);
            
            % Compute the shear force components
            qShear = BShearForces*dHatActual;
            
            % Decide upon the the visualization resultant
            if strcmp(graph.component,'1')
                resultantComponent(xiCounter,etaCounter) = qShear(1,1);
            elseif strcmp(graph.component,'2')
                resultantComponent(xiCounter,etaCounter) = qShear(2,1);
            end
        end

        % Update the counter in xi-direction
        xiCounter = xiCounter + 1;
        
        % Update the parametric coordinate in xi-direction
        xi = xi + dxi;
    end
    % Update the counter in eta-direction
    etaCounter = etaCounter + 1;
    
    % Update the parametric coordinate in eta-direction
    eta = eta + deta;
end

%% 2. Visualize the resultant and the knots over the domain

% Plot the resultant over the reference configuration
surf(P(:,:,1),P(:,:,2),P(:,:,3),resultantComponent(:,:));
hold on;

% Plot the element boundaries on the undeformed geometry
plot_knotsForBSplineSurfaceOnCartesianSpace(p,q,Xi,Eta,CP,isNURBS,isUndeformed,xiGrid,etaGrid);
hold off;

end
