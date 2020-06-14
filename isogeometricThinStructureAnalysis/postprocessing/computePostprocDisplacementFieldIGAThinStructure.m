function d = computePostprocDisplacementFieldIGAThinStructure ...
    (BSplinePatch, dHat, xi, eta)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the displacement field at a given parametric location of a
% BSpline patch given the Control Point displacement field of the patch.
%
%        Input :
% BSplinePatch : Parameters of the B-Spline patch :
%                   .p,.q : Polynomial orders along xi and eta parametric 
%                           directions
%                .Xi,.Eta : Knot vectors along xi and eta parametric 
%                           directions
%                     .CP : Set of Control Point coordinates and weights
%                .isNURBS : Flag on whether the basis functions are
%                           B-Splines or NURBS
%         dHat : The Control Point displacement vector
%       xi,eta : The parametric locations where to compute the displacement
%                field
%               
%       Output :
%            d : The displacement field at the given parametric location
%                d = d(xi,eta)
%
% Function main body :
%
% 0. Read input
%
% 1. Split the control point displacement field to each element
%
% 2. Compute the basis functions at the given parametric location
%
% 3. Compute the control point displacement field at the element where the given parametric coordinates belong to
%
% 4. Compute the displacement field at the given parametric location
%
%% Function main body

%% 0. Read input
p = BSplinePatch.p;
Xi = BSplinePatch.Xi;
q = BSplinePatch.q;
Eta = BSplinePatch.Eta;
CP = BSplinePatch.CP;
isNURBS = BSplinePatch.isNURBS;
mxi = length(Xi);
meta = length(Eta);
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));

%% 1. Split the control point displacement field to each element
dHatElmnt = zeros(mxi - p - 1,meta - q - 1,3*(p + 1)*(q + 1));
for iEtaSpan = (q + 1):(meta - q - 1)
    for iXiSpan = (p + 1):(mxi - p - 1)
        xiCounter = 1; 
        for c = iEtaSpan - q - 1:iEtaSpan - 1 
            for b = iXiSpan - p:iXiSpan
                dHatElmnt(iXiSpan,iEtaSpan,xiCounter) = dHat(3*(c*nxi + b) - 2);
                dHatElmnt(iXiSpan,iEtaSpan,xiCounter + 1) = dHat(3*(c*nxi + b) - 1);
                dHatElmnt(iXiSpan,iEtaSpan,xiCounter + 2) = dHat(3*(c*nxi + b));
                xiCounter = xiCounter + 3;
            end
        end
    end
end

%% 2. Compute the basis functions at the given parametric location
xiSpan = findKnotSpan(xi,Xi,nxi);
etaSpan = findKnotSpan(eta,Eta,neta);
R = computeIGABasisFunctionsAndDerivativesForSurface...
    (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,0);

%% 3. Compute the control point displacement field at the element where the given parametric coordinates belong to
dHatActual = squeeze(dHatElmnt(xiSpan,etaSpan,:));

%% 4. Compute the displacement field at the given parametric location
d = computePostprocDisplacementIGAKirchhoffLoveShell...
    (p,q,R,dHatActual);

end
