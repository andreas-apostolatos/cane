function homDOFs = findDOFsClampedSupportIGAKLShell...
    (homDOFs, xiSup, etaSup, p, Xi, q, Eta, CP, isNURBS)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the global numbering of the DOFs which need to be contrained such
% that to obtain a clamped support at one of the 4 edges of the
% isogeometric Kirchhoff-Love shell.
%
%       Input :
%     homDOFs : The outdated vector of the DOFs where homogeneous Dirichlet
%               boundary conditions are applied
%       xiSup : The support parametric extension in xi-direction
%      etaSup : The support parametric extension in eta-direction
%           p : The polynomial order in xi-direction
%          Xi : The knot vector in xi-directions
%           q : The polynomial order in eta-direction
%         Eta : The knot vector in eta-directions
%          CP : The set of Control Points and weights
%     isNURBS : Flag on whether the basis functions are B-Spline or NURBS
%
%      Output :
%     homDOFs : The updated vector of the DOFs where homogeneous Dirichlet
%               boundary conditions are applied
%
% Function layout :
%
% 0. Read input
%
% 1. Find the middle point of the edge
%
% 2. Compute the normal to the surface basis vector at the middle point of the edge
%
% 3. Determine the orientation of the clamping
%
% 4. Get the global numbering of the DOFs which need to be fixed for a fixed support
%
% 5. Get the global numbering of the DOFs which need to be fixed in order to extend a fixed to a clamped support
%
% 6. Sort the DOFs and delete dublicate values
%
%% Function main body

%% 0. Read input

% Assign a tolerance value
eps = 1e-12;

% Number of knots in xi-direction
mxi = length(Xi);

% Number of Control Points in xi-direction
nxi = length(CP(:,1,1));

% Number of knots in eta-direction
meta = length(Eta);

% Number of Control Points in eta-direction
neta = length(CP(1,:,1));

% Check input
checkInputForBSplineSurface(p,mxi,nxi,q,meta,neta);

%% 1. Find the middle point of the edge
if xiSup(1) == xiSup(2)
    xi = xiSup(1);
    eta = (etaSup(1) + etaSup(2))/2;
    isOnXi = false;
elseif etaSup(1) == etaSup(2)
    xi = (xiSup(1) + xiSup(2))/2;
    eta = etaSup(1);
    isOnXi = true;
else
    error('The given parametric extensions do not correspond to a boundary of the patch');
end

%% 2. Compute the normal to the surface basis vector at the middle point of the edge

% Find the knot span indices
xiSpan = findKnotSpan(xi,Xi,nxi);
etaSpan = findKnotSpan(eta,Eta,neta);

% Compute the basis functions and their derivatives
dR = computeIGABasisFunctionsAndDerivativesForSurface...
    (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,1);

% Compute the base vectors of the surface
[A1,A2] = computeBaseVectorsAndDerivativesForBSplineSurface...
    (xiSpan,p,etaSpan,q,CP,0,dR);

% Compute the surface normal base vector
A3 = cross(A1,A2);
A3 = A3/norm(A3);

%% 3. Determine the orientation of the clamping
if abs(1 - abs(A3(1,1))) <= eps && abs(A3(2,1)) <= eps && abs(A3(3,1)) <= eps
    dirSup = 1;
elseif abs(A3(1,1)) <= eps && abs(1 - abs(A3(2,1))) <= eps && abs(A3(3,1)) <= eps
    dirSup = 2;
elseif abs(A3(1,1)) <= eps && abs(A3(2,1)) <= eps && abs(1 - abs(A3(3,1))) <= eps
    dirSup = 3;
else
    error('The orientation of the surface at the clamping support is not aligned with the Cartesian axes');
end

%% 4. Get the global numbering of the DOFs which need to be fixed for a fixed support
homDOFs = findDofs3D(homDOFs,xiSup,etaSup,1,CP);
homDOFs = findDofs3D(homDOFs,xiSup,etaSup,2,CP);
homDOFs = findDofs3D(homDOFs,xiSup,etaSup,3,CP);

%% 5. Get the global numbering of the DOFs which need to be fixed in order to extend a fixed to a clamped support
noFixedDOFs = length(homDOFs);
if ~isOnXi && xi == 0
    for j = 1:neta
        homDOFs(noFixedDOFs + j) = (3 + dirSup) + 3*nxi*(j-1);
    end
elseif isOnXi && eta == 0
    for j = 1:nxi
        homDOFs(noFixedDOFs + j) = 3*(nxi + j) - (3 - dirSup);
    end
elseif isOnXi && eta == 1
    for j = 1:nxi
        homDOFs(noFixedDOFs + j) = 3*(nxi*(neta-2) + j) - (3 - dirSup);
    end
end

%% 6. Sort the DOFs and delete dublicate values
homDOFs = sort(homDOFs);
homDOFs = unique(homDOFs);

end
