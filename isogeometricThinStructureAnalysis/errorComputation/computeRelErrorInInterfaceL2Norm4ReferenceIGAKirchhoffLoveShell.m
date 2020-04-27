function [displRelDiffL2, forceRelDiffL2, momentRelDiffL2] = ...
    computeRelErrorInInterfaceL2Norm4ReferenceIGAKirchhoffLoveShell ...
    (patch, dHat, patchReference, dHatReference, newtonRaphson, intC)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the relative error in the L2-, H1-norm for the displacement field
% , the L2 norm for the rotational field, the force tensor field, the 
% moment tensor field and the shear force traction field over the domain 
% for an isogeometric Kirchhoff-Love shell patch compared to a reference 
% solution which is also provided in NURBS parametrization. The L2-norm is 
% computed over the boundary of the patch which is assumed to be one of 
% those which are coupled to solve the domain decomposed elasticity 
% boundary value problem.
%
%           Input :
%           patch : The isogeometric Kirchhoff-Love shell patch for which
%                   the relative error to be computed
%            dHat : The numerical solution over the given Kirchhoff-Love 
%                   shell patch
%  patchReference : The single patch model over which the reference solution
%                   is computed
%   dHatReference : The reference solution over the single patch model
%   newtonRaphson :   newtonRaphson.eps : Residual tolerance
%                   newtonRaphson.maxIt : Maximum number of iterations
%            intC : On the numerical integration over line
%
%          Output :
%  displRelDiffL2 : The relative error of the displacement field in the
%                   L2-norm
%  forceRelDiffL2 : The relative error of the force tensor field in the
%                   L2-norm
% momentRelDiffL2 : The relative error of the momnet tensor field in the
%                   L2-norm
%
% Function layout :
%
% 0. Read input
%
% 1. Distribute the nodal displacements into the elements
%
% 2. Loop over all interfaces involving the given patch
% ->
%    2i. Get the interface parametrization
%
%   2ii. Get the running and the fixed parameters on the patch interface and the integration region over the knot vector
%
%  2iii. Issue Gauss Point coordinates and weights
%
%   2iv. Loop over all the elements on the integration interface
%   ->
%        2iv.1. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
%
%        2iv.2. Loop over all Gauss points
%        ->
%               2iv.1. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
%
%              2iv.2ii. Compute the NURBS basis functions for the patch
%
%             2iv.2iii. Compute the physical location of the Gauss Point on patch
%
%              2iv.2iv. Compute the parameter location of the Gauss Point on patch in patchReference
%
%               2iv.2v. Find the knot span indices for the reference patch
%
%              2iv.2vi. Compute the NURBS basis functions for the reference patch
%
%             2iv.2vii. Get the actual Control Point displacement vectors for the current knot span
%               
%            2iv.2viii. Compute the base vectors of the configuration over the interface for the patch
%
%              2iv.2ix. Compute the actual displacement vectors at the Gauss point
%
%               2iv.2x. Compute the membrane force tensor at the Gauss point
%
%              2iv.2xi. Compute the moment tensor at the Gauss point
%
%             2iv.2xii. Compute the normal and the tangent to the boundary vector
%
%            2iv.2xiii. Compute the element length at the GP
%
%             2iv.2xiv. Compute the L2-norms
%        <-
%   <-
% <-
% 3. Get the relative values
%
%% Function main body

%% 0. Read input

% For patch 1 :
% _____________

% Reassign the analysis arrays
p = patch.p;
q = patch.q;
Xi = patch.Xi;
Eta = patch.Eta;
CP = patch.CP;
isNURBS = patch.isNURBS;
xicoup = patch.xicoup;
etacoup = patch.etacoup;
parameters = patch.parameters;

% Check number of interfaces involving the given patch
noInterfacesXi = length(xicoup(:,1));
noInterfacesEta = length(etacoup(:,1));
if noInterfacesXi ~= noInterfacesEta
    error('No. of interfaces in xi and in eta direction does not match');
else
    noInterfaces = noInterfacesXi;
end

% Number of knots in u,v-direction
mxi = length(Xi);
meta = length(Eta);

% Number of Control Points in xi-,eta- directions
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

% For patch 2 :
% _____________

pRef = patchReference.p;
qRef = patchReference.q;
XiRef = patchReference.Xi;
EtaRef = patchReference.Eta;
CPRef = patchReference.CP;
isNURBSRef = patchReference.isNURBS;
parametersRef = patchReference.parameters;

% Number of knots in u,v-direction
mxiRef = length(XiRef);
metaRef = length(EtaRef);

% Number of Control Points in xi-,eta- directions
nxiRef = length(CPRef(:,1,1));
netaRef = length(CPRef(1,:,1));

% Compute the membrane material matrix
DmRef = parametersRef.E*parametersRef.t/(1-parametersRef.nue^2)*...
        [1              parametersRef.nue 0
         parametersRef.nue 1              0
         0              0                 (1-parametersRef.nue)/2];
                                                 
% Compute the bending material matrix
DbRef = parametersRef.E*parametersRef.t^3/(12*(1-parametersRef.nue^2))*...
        [1              parametersRef.nue 0
         parametersRef.nue 1              0 
         0              0                 (1-parametersRef.nue)/2];
     
% Define a tolerance
tolerance = 1e-12;
     
% On the Netwon-Raphson iteration parameters
xi0 = .5;
eta0 = .5;

% Initialize auxiliary arrays
forceMatrix = zeros(2,2);
forceMatrixRef = zeros(2,2);
momentMatrix = zeros(2,2);
momentMatrixRef = zeros(2,2);

% Initialize auxiliary variables
dKommaX = zeros(3,1);
dKommaY = zeros(3,1);
dKommaZ = zeros(3,1);

% Initialize output values
displDiffL2 = 0;
displRefL2 = 0;
forceDiffL2 = 0;
forceRefL2 = 0;
momentDiffL2 = 0;
momentRefL2 = 0;

%% 1. Distribute the nodal displacements into the elements

% For patch :
% ___________

dHatElem = zeros(mxi-p-1,meta-q-1,3*(p+1)*(q+1));
for cEtaSpan = (q+1):(meta-q-1)
    for cXiSpan = (p+1):(mxi-p-1)
        xiCounter = 1; 
        for c = cEtaSpan-q-1:cEtaSpan-1 
            for b = cXiSpan-p:cXiSpan
                dHatElem(cXiSpan,cEtaSpan,xiCounter) = dHat(3*(c*nxi+b)-2);
                dHatElem(cXiSpan,cEtaSpan,xiCounter + 1) = dHat(3*(c*nxi+b)-1);
                dHatElem(cXiSpan,cEtaSpan,xiCounter + 2) = dHat(3*(c*nxi+b));

                % Update counter
                xiCounter = xiCounter + 3;
            end
        end
    end
end

% For patchReference :
% ____________________

dHatElemRef = zeros(mxiRef-pRef-1,metaRef-qRef-1,3*(pRef+1)*(qRef+1));
for cEtaSpan = (qRef+1):(metaRef-qRef-1)
    for cXiSpan = (pRef+1):(mxiRef-pRef-1)
        xiCounter = 1; 
        for c = cEtaSpan-qRef-1:cEtaSpan-1 
            for b = cXiSpan-pRef:cXiSpan
                dHatElemRef(cXiSpan,cEtaSpan,xiCounter) = dHatReference(3*(c*nxiRef+b)-2);
                dHatElemRef(cXiSpan,cEtaSpan,xiCounter + 1) = dHatReference(3*(c*nxiRef+b)-1);
                dHatElemRef(cXiSpan,cEtaSpan,xiCounter + 2) = dHatReference(3*(c*nxiRef+b));

                % Update counter
                xiCounter = xiCounter + 3;
            end
        end
    end
end

%% 2. Loop over all interfaces involving the given patch
for cInter = 1:noInterfaces
    %% 2i. Get the interface parametrization
    xicoupPatch = xicoup(cInter,:);
    etacoupPatch = etacoup(cInter,:);
    
    %% 2ii. Get the running and the fixed parameters on the patch interface and the integration region over the knot vector
    if etacoupPatch(1) == etacoupPatch(2)
        % Coupled region in xi-direction
        intRegion = xicoupPatch;

        % Find the correct spans for the coupled region
        spanStart = findKnotSpan(intRegion(1),Xi,nxi);
        spanEnd = findKnotSpan(intRegion(2),Xi,nxi) + 1;

        % Corresponding to the coupled region knot span
        intRegionOnKnotVector = Xi(spanStart:spanEnd);

        % Fixed parameter on the parametric net
        eta = etacoupPatch(1);

        % Find the span where uvL it lies in
        etaSpan = findKnotSpan(eta,Eta,neta);

        % Flag on whether the coupling line is over xi
        isOnXi = true;
        
        % On the initial guess for the Newton Raphson iterations
%         eta0 = eta;
    elseif xicoupPatch(1) == xicoupPatch(2)
        % Coupled region in eta-direction
        intRegion = etacoupPatch;

        % Find the correct spans for the coupled region
        spanStart = findKnotSpan(intRegion(1),Eta,neta);   
        spanEnd = findKnotSpan(intRegion(2),Eta,neta) + 1;   

        % Corresponding to the coupled region knot span
        intRegionOnKnotVector = Eta(spanStart:spanEnd);

        % Fixed parameter on the parametric net
        xi = xicoupPatch(1);

        % Find the span where uv it lies in
        xiSpan = findKnotSpan(xi,Xi,nxi);

        % Flag on whether the coupling line is over eta
        isOnXi = false;
        
        % On the initial guess for the Newton Raphson iterations
%         xi0 = xi;
    end

    %% 2iii. Issue Gauss Point coordinates and weights
    if strcmp(intC.type,'default')
        if isOnXi
            nGPs = p + 1;
        else
            nGPs = q + 1;
        end
    else
        nGPs = intC.nGPError;
    end
    [GP,GW] = getGaussPointsAndWeightsOverUnitDomain(nGPs);

    %% 2iv. Loop over all the elements on the integration interface
    for i=1:length(intRegionOnKnotVector)-1
        if intRegionOnKnotVector(i) ~= intRegionOnKnotVector(i+1)
            %% 2iv.1. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
            detJxizeta = (intRegionOnKnotVector(i+1)-intRegionOnKnotVector(i))/2;

            %% 2iv.2. Loop over all Gauss points
            for j=1:nGPs
                %% 2iv.2i. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
                xiEta = ((1-GP(j))*intRegionOnKnotVector(i)+(1+GP(j))*intRegionOnKnotVector(i+1))/2;

                %% 2iv.2ii. Compute the NURBS basis functions for the patch
                if isOnXi
                    xi = xiEta;
                    xiSpan = findKnotSpan(xi,Xi,nxi);
                else
                    eta = xiEta;
                    etaSpan = findKnotSpan(eta,Eta,neta);
                end
                dR = computeIGABasisFunctionsAndDerivativesForSurface(xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,2);

                %% 2iv.2iii. Compute the physical location of the Gauss Point on patch
                P = computeCartesianCoordinatesOfAPointOnBSplineSurface(xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,dR(:,1));

                %% 2iv.2iv. Compute the parameter location of the Gauss Point on patch in patchReference
                [xiRef,etaRef,PTest,flagNR,~] = computeNearestPointProjectionOnBSplineSurface...
                    (P,pRef,XiRef,qRef,EtaRef,CPRef,isNURBSRef,xi0,eta0,newtonRaphson);
                if ~flagNR
                    error('Closest point projection on NURBS surface failed');
                end
                
                %% 2iv.2v. Find the knot span indices for the reference patch
                xiSpanRef = findKnotSpan(xiRef,XiRef,nxiRef);
                etaSpanRef = findKnotSpan(etaRef,EtaRef,netaRef);
                
                %% 2iv.2vi. Compute the NURBS basis functions for the reference patch
                dRRef = computeIGABasisFunctionsAndDerivativesForSurface(xiSpanRef,pRef,xiRef,XiRef,etaSpanRef,qRef,etaRef,EtaRef,CPRef,isNURBSRef,2);
            
                %% 2iv.2vii. Get the actual Control Point displacement vectors for the current knot span

                % For patch :
                % ___________

                dHatActual(:,1) = dHatElem(xiSpan,etaSpan,:);

                % For patchReference :
                % ____________________

                dHatActualRef(:,1) = dHatElemRef(xiSpanRef,etaSpanRef,:);

                %% 2iv.2viii. Compute the base vectors of the configuration over the interface for the patch

                % For patch :
                % ___________

                % Compute the in-plane base vectors
                [dG1,dG2] = ...
                    computeBaseVectorsAndDerivativesForBSplineSurface...
                    (xiSpan,p,etaSpan,q,CP,1,dR);

                % Compute the not normalized surface normal vector
                G3Tilde = cross(dG1(:,1),dG2(:,1));

                % Compute the Jacobian of the transformation from the physical
                % to the parameter space
                detJPhys2Param = norm(G3Tilde);

                % Compute the surface normal base vector
                G3 = G3Tilde/detJPhys2Param;

                % For patchReference :
                % ____________________

                % Compute the in-plane base vectors
                [dG1Ref,dG2Ref] = ...
                    computeBaseVectorsAndDerivativesForBSplineSurface...
                    (xiSpanRef,pRef,etaSpanRef,qRef,CPRef,1,dRRef);

                % Compute the not normalized surface normal vector
                G3TildeRef = cross(dG1Ref(:,1),dG2Ref(:,1));

                % Compute the Jacobian of the transformation from the physical
                % to the parameter space
                detJPhys2ParamRef = norm(G3TildeRef);

                % Compute the surface normal base vector
                G3Ref = G3TildeRef/detJPhys2ParamRef;

                %% 2iv.2ix. Compute the actual displacement vectors at the Gauss point

                % For patch :
                % ___________

                d = computePostprocDisplacementIGAKirchhoffLoveShell(p,q,dR(:,1),dHatActual);

                % For patchReference :
                % ____________________

                dRef = computePostprocDisplacementIGAKirchhoffLoveShell(pRef,qRef,dRRef(:,1),dHatActualRef);

                %% 2iv.2x. Compute the membrane force tensor at the Gauss point

                % For patch :
                % ___________

                % Compute the B-operator matrix for the strain field
                BStrain = computeBOperatorMatrix4StrainIGAKirchhoffLoveShellLinear...
                    (p,q,dR,[dG1(:,1) dG2(:,1)]);

                % Compute the strain tensor in a Voigt notation
                epsilonVoigt = BStrain*dHatActual;

                % Compute the force tensor in a Voigt notation
                forceVoigt = Dm*epsilonVoigt;

                % Assemble the force tensor in a matrix notation
                forceMatrix(1,1) = forceVoigt(1,1);
                forceMatrix(1,2) = forceVoigt(3,1);
                forceMatrix(2,1) = forceMatrix(1,2);
                forceMatrix(2,2) = forceVoigt(2,1);

                % For patchReference :
                % ____________________

                % Compute the B-operator matrix for the strain field
                BStrainRef = computeBOperatorMatrix4StrainIGAKirchhoffLoveShellLinear...
                    (pRef,qRef,dRRef,[dG1Ref(:,1) dG2Ref(:,1)]);

                % Compute the strain tensor in a Voigt notation
                epsilonVoigtRef = BStrainRef*dHatActualRef;

                % Compute the force tensor in a Voigt notation
                forceVoigtRef = DmRef*epsilonVoigtRef;

                % Assemble the force tensor in a matrix notation
                forceMatrixRef(1,1) = forceVoigtRef(1,1);
                forceMatrixRef(1,2) = forceVoigtRef(3,1);
                forceMatrixRef(2,1) = forceMatrixRef(1,2);
                forceMatrixRef(2,2) = forceVoigtRef(2,1);

                %% 2iv.2xi. Compute the moment tensor at the Gauss point

                % For patch :
                % ___________

                % Compute the B-operator matrix for the change in the curvature 
                % field
                BCurvature = computeBOperatorMatrix4CurvatureIGAKirchhoffLoveShellLinear(p,q,dR,[dG1(:,1) dG2(:,1)],[dG1(:,2) dG2(:,2) dG1(:,3)],G3Tilde);

                % Compute the strain tensor in a Voigt notation
                curvatureVoigt = BCurvature*dHatActual;

                % Compute the force tensor in a Voigt notation
                momentVoigt = Db*curvatureVoigt;

                % Assemble the force tensor in a matrix notation
                momentMatrix(1,1) = momentVoigt(1,1);
                momentMatrix(1,2) = momentVoigt(3,1);
                momentMatrix(2,1) = momentMatrix(1,2);
                momentMatrix(2,2) = momentVoigt(2,1);

                % For patchReference :
                % ____________________

                % Compute the B-operator matrix for the change in the curvature 
                % field
                BCurvatureRef = computeBOperatorMatrix4CurvatureIGAKirchhoffLoveShellLinear(pRef,qRef,dRRef,[dG1Ref(:,1) dG2Ref(:,1)],[dG1Ref(:,2) dG2Ref(:,2) dG1Ref(:,3)],G3TildeRef);

                % Compute the strain tensor in a Voigt notation
                curvatureVoigtRef = BCurvatureRef*dHatActualRef;

                % Compute the force tensor in a Voigt notation
                momentVoigtRef = DbRef*curvatureVoigtRef;

                % Assemble the force tensor in a matrix notation
                momentMatrixRef(1,1) = momentVoigtRef(1,1);
                momentMatrixRef(1,2) = momentVoigtRef(3,1);
                momentMatrixRef(2,1) = momentMatrixRef(1,2);
                momentMatrixRef(2,2) = momentVoigtRef(2,1);


                %% 2iv.2xii. Compute the normal and the tangent to the boundary vector
                [n,t] = computeNormalAndTangentVectorsToBSplineBoundary(xi,Xi,eta,Eta,dG1(:,1),dG2(:,1),G3,isOnXi);

                %% 2iv.2xiii. Compute the element length at the GP
                elementLengthOnGP = detJPhys2Param*detJxizeta*GW(j);

                %% 2iv.2xiv. Compute the L2-norms

                % Compute the L2-norms of the differences

                % Displacement vector
                displDiffL2 = displDiffL2 + norm(d - dRef)^2*elementLengthOnGP;

                % Force tensor
                forceDiffL2 = forceDiffL2 + norm(forceMatrix - forceMatrixRef)^2*elementLengthOnGP;

                % Moment tensor
                momentDiffL2 = momentDiffL2 + norm(momentMatrix - momentMatrixRef)^2*elementLengthOnGP;

                % Compute the L2-norms of the reference fields

                % Displacement
                displRefL2 = displRefL2 + norm(dRef)^2*elementLengthOnGP;

                % Force tensor
                forceRefL2 = forceRefL2 + norm(forceMatrixRef)^2*elementLengthOnGP;

                % Moment tensor
                momentRefL2 = momentRefL2 + norm(momentMatrixRef)^2*elementLengthOnGP;
            end
        end
    end
end

%% 3. Get the relative values

% Displacement vector
if sqrt(displRefL2) > tolerance
    displRelDiffL2 = sqrt(displDiffL2)/sqrt(displRefL2);
else
    displRelDiffL2 = sqrt(displDiffL2);
end

% Force tensor
if sqrt(forceRefL2) > tolerance
    forceRelDiffL2 = sqrt(forceDiffL2)/sqrt(forceRefL2);
else
    forceRelDiffL2 = sqrt(forceDiffL2);
end

% Moment tensor
if sqrt(momentRefL2) > tolerance
    momentRelDiffL2 = sqrt(momentDiffL2)/sqrt(momentRefL2);
else
    momentRelDiffL2 = sqrt(momentDiffL2);
end

end
