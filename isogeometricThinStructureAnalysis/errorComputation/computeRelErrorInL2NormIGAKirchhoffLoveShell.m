function [displRelDiffL2, forceRelDiffL2, momentRelDiffL2] = ...
    computeRelErrorInL2NormIGAKirchhoffLoveShell ...
    (patch, dHat, patchReference, dHatReference, newtonRaphson, int)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Reuturns the error in the L2-norm for the displacement, the force and the
% moment tensor in the shell's domain, given the solution obtained from a
% refined patch taken as a reference solution in terms of displacement for
% the Kirchhoff-Love shell problem.
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
%             int : On the numerical integration over line
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
% 1. Compute the integration knot spans
%
% 2. Distribute the nodal displacements into the elements
%
% 3. Choose the Gauss Integration rule
%
% 4. Loop over the elements
% ->
%    4i. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
%
%   4ii. Loop over the Gauss points
%   ->
%        4ii.1. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
%
%        4ii.2. Find the knot span
%
%        4ii.3. Compute the NURBS basis functions for the patch
%
%        4ii.4. Compute the physical location of the Gauss Point on patch
%
%        4ii.5. Compute the parameter location of the Gauss Point on patch in patchReference
%
%        4ii.6. Find the knot span indices for the reference patch
%
%        4ii.7. Compute the NURBS basis functions for the reference patch
%
%        4ii.8. Get the actual Control Point displacement vectors for the current knot span
%
%        4ii.9. Compute the base vectors of the configuration over the interface for the patch
%
%        4ii.10. Compute the actual displacement vectors at the Gauss point
%
%        4ii.11. Compute the force tensor at the Gauss point
%
%        4ii.12. Compute the moment tensor at the Gauss point
%
%        4ii.13. Compute the element length at the GP
%
%        4ii.14. Compute the L2-norms
%   <-
% <-
%
% 5. Get the relative values
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
parameters = patch.parameters;

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

% For the reference patch :
% _________________________

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
eps1 = 1e-9;
eps2 = 1e-9;
maxIt = 5;

% Initialize auxiliary arrays
forceTensor = zeros(2,2);
forceTensorRef = zeros(2,2);
momentTensor = zeros(2,2);
momentTensorRef = zeros(2,2);

% Initialize auxiliary variables
dKommaX = zeros(3,1);
dKommaY = zeros(3,1);
dKommaZ = zeros(3,1);

% Initialize the reference values for the displacement vector, the force 
% and the moment tensor if those are not provided

% Initialize output values
displDiffL2 = 0;
displRefL2 = 0;
forceDiffL2 = 0;
forceRefL2 = 0;
momentDiffL2 = 0;
momentRefL2 = 0;

%% 1. Compute the integration knot spans

% Integration knot span in xi-direction
XiMerged = mergesorted(Xi,XiRef);
XiMerged = unique(XiMerged);

% Integration knot span in eta-direction
EtaMerged = mergesorted(Eta,EtaRef);
EtaMerged = unique(EtaMerged);

%% 2. Distribute the nodal displacements into the elements

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

%% 3. Choose the Gauss Integration rule
if strcmp(int.type,'default')
    noGPsXi = p + 1;
    noGPsEta = q + 1;
elseif strcmp(int.type,'user')
    noGPsXi = int.noGPXiError;
    noGPsEta = int.noGPEtaError;
end
[GPXi,GWXi] = getGaussPointsAndWeightsOverUnitDomain(noGPsXi);
[GPEta,GWEta] = getGaussPointsAndWeightsOverUnitDomain(noGPsEta);

%% 4. Loop over the elements
for counterEta = 1:length(unique(Eta)) - 1
    for counterXi = 1:length(unique(Xi)) - 1
        %% 4i. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
        detJxizeta = (XiMerged(counterXi+1)-XiMerged(counterXi))*(EtaMerged(counterEta+1)-EtaMerged(counterEta))/4;
                
        %% 4ii. Loop over the Gauss points
        for counterGPEta = 1:length(GWEta)
            for counterGPXi = 1:length(GWXi)
                %% 4ii.1. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
                xi = ((1-GPXi(counterGPXi))*XiMerged(counterGPXi) + (1+GPXi(counterGPXi))*XiMerged(counterGPXi+1))/2;
                eta = ((1-GPEta(counterGPEta))*EtaMerged(counterGPEta) + (1+GPEta(counterGPEta))*EtaMerged(counterGPEta+1))/2;
                
                %% 4ii.2. Find the knot span
                xiSpan = findKnotSpan(xi,Xi,nxi);
                etaSpan = findKnotSpan(eta,Eta,neta);
                
                %% 4ii.3. Compute the NURBS basis functions for the patch
                noDeriv = 2;
                dR = computeIGABasisFunctionsAndDerivativesForSurface...
                    (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,noDeriv);

                %% 4ii.4. Compute the physical location of the Gauss Point on patch
                P = computeCartesianCoordinatesOfAPointOnBSplineSurface...
                    (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,dR(:,1));

                %% 4ii.5. Compute the parameter location of the Gauss Point on patch in patchReference
                [xiRef,etaRef,~,flagNR,~] = computeNearestPointProjectionOnBSplineSurface...
                    (P,pRef,XiRef,qRef,EtaRef,CPRef,isNURBSRef,xi0,eta0,newtonRaphson);
                if ~flagNR
                    error('Closest point projection on NURBS surface failed');
                end

                %% 4ii.6. Find the knot span indices for the reference patch
                xiSpanRef = findKnotSpan(xiRef,XiRef,nxiRef);
                etaSpanRef = findKnotSpan(etaRef,EtaRef,netaRef);

                %% 4ii.7. Compute the NURBS basis functions for the reference patch
                noDeriv = 2;
                dRRef = computeIGABasisFunctionsAndDerivativesForSurface...
                    (xiSpanRef,pRef,xiRef,XiRef,etaSpanRef,qRef,etaRef,EtaRef,CPRef,isNURBSRef,noDeriv);

                %% 4ii.8. Get the actual Control Point displacement vectors for the current knot span

                % For patch :
                % ___________

                dHatActual(:,1) = dHatElem(xiSpan,etaSpan,:);

                % For patchReference :
                % ____________________

                dHatActualRef(:,1) = dHatElemRef(xiSpanRef,etaSpanRef,:);

                %% 4ii.9. Compute the base vectors of the configuration over the interface for the patch

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

                %% 4ii.10. Compute the actual displacement vectors at the Gauss point

                % For patch :
                % ___________

                d = computePostprocDisplacementIGAKirchhoffLoveShell(p,q,dR(:,1),dHatActual);

                % For patchReference :
                % ____________________

                dRef = computePostprocDisplacementIGAKirchhoffLoveShell...
                    (pRef,qRef,dRRef(:,1),dHatActualRef);

                %% 4ii.11. Compute the force tensor at the Gauss point

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
                forceTensor(1,1) = forceVoigt(1,1);
                forceTensor(1,2) = forceVoigt(3,1);
                forceTensor(2,1) = forceTensor(1,2);
                forceTensor(2,2) = forceVoigt(2,1);

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
                forceTensorRef(1,1) = forceVoigtRef(1,1);
                forceTensorRef(1,2) = forceVoigtRef(3,1);
                forceTensorRef(2,1) = forceTensorRef(1,2);
                forceTensorRef(2,2) = forceVoigtRef(2,1);

                %% 4ii.12. Compute the moment tensor at the Gauss point

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
                momentTensor(1,1) = momentVoigt(1,1);
                momentTensor(1,2) = momentVoigt(3,1);
                momentTensor(2,1) = momentTensor(1,2);
                momentTensor(2,2) = momentVoigt(2,1);

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
                momentTensorRef(1,1) = momentVoigtRef(1,1);
                momentTensorRef(1,2) = momentVoigtRef(3,1);
                momentTensorRef(2,1) = momentTensorRef(1,2);
                momentTensorRef(2,2) = momentVoigtRef(2,1);

                %% 4ii.13. Compute the element length at the GP
                elementAreaOnGP = detJPhys2Param*detJxizeta*GWXi(counterGPXi)*GWXi(counterGPEta);

                %% 4ii.14. Compute the L2-norms

                % Compute the L2-norms of the differences

                % Displacement vector
                displDiffL2 = displDiffL2 + norm(d - dRef)^2*elementAreaOnGP;

                % Force tensor
                forceDiffL2 = forceDiffL2 + norm(forceTensor - forceTensorRef)^2*elementAreaOnGP;

                % Moment tensor
                momentDiffL2 = momentDiffL2 + norm(momentTensor - momentTensorRef)^2*elementAreaOnGP;

                % Compute the L2-norms of the reference fields

                % Displacement
                displRefL2 = displRefL2 + norm(dRef)^2*elementAreaOnGP;

                % Force tensor
                forceRefL2 = forceRefL2 + norm(forceTensorRef)^2*elementAreaOnGP;

                % Moment tensor
                momentRefL2 = momentRefL2 + norm(momentTensorRef)^2*elementAreaOnGP;
            end
        end
    end
end

%% 5. Get the relative values

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
