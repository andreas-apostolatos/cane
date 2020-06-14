function [F, tanMtx] = computeLoadVctPointIGAThinStructure ...
    (FOutdated, BSplinePatch, xib, etab, FAmp, direction, isFollower, t, ...
    propInt, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the consistent force vector for a surface B-Spline patch
% corresponding to a surface load.
%
%          Input :
%      FOutdated : The existing force vector
%   BSplinePatch : Array containing information on the B-Spline patch
%                     .p,.q : Polynomial orders in xi- and -eta direction
%                  .Xi,.Eta : Polynomial orders in xi-eta direction
%                       .CP : Set of Control Point coordinates and weights
%                  .isNURBS : Flag on whether the basis is a NURBS or a
%                             B-Spline
%                   .noDOFs : Total Number of DOFs for the B-Spline patch
%                             including possible Lagrange Multipliers
%       xib,etab : Load position in the B-Spline parameter space
%           FAmp : The magnitude of the point load
%      direction : Direction of the applied force
%                           'x' : along x -> FAmp*e1 
%                           'y' : along y -> FAmp*e2
%                           'z' : along z -> FAmp*e3
%                      'theta1' : along theta1 -> FAmp*G1/norm(G1)
%                      'theta2' : along theta2 -> FAmp*G2/norm(G2)
%                      'normal' : normal to the surface -> FAmp*n
%     isFollower : Flag on whether the applied load is a follower load
%        propInt : On the integral integration (dummy variable for this 
%                  function)
%         outMsg : Whether or not to output message
%                  'outputEnabled' : enables output information   
%
%         Output :
%             Fl : The updated force vector
%         tanMtx : The tangent matrix resulting from the application of a
%                  follower load
%
% Function layout :
%
% 0. Read input
%
% 1. Find the knot span of the parametric location where the load is applied
%
% 2. Get the Element Freedom Table (EFT) corresponding to the element where the load application parametric location belongs
%
% 3. Compute the IGA basis functions and their derivatives at the parametric location where the load is applied
%
% 4. Compute the base vectors
%
% 5. Compute the normal to the surface vector
%
% 6. Compute the traction vector
%
% 7. Compute the matric containing the basis functions
%
% 8. Compute the consistent load vector corresponding to the point load application
%
% 9. Update by the existing load vector
%
% 10. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('_______________________________________________________\n');
    fprintf('#######################################################\n');
    if isvector(FOutdated)
        fprintf('Update of the load vector corresponding to point\n');
    else
        fprintf('Computation of the load vector corresponding to point\n');
    end
    fprintf('load for the isogeometric Kirchhoff-Love shell problem\n');
    fprintf('has been initiated\n\n');
    if strcmp(direction,'x') || strcmp(direction,'y') || strcmp(direction,'z') ...
            || strcmp(direction,'theta1') || strcmp(direction,'theta2') ...
            || strcmp(direction,'normal')
        fprintf('Load direction is chosen as %s',direction);
    else
        error('Select a valid direction for the load');
    end
    fprintf('\n');
    if isFollower
        fprintf('Follower load is considered');
    else
        fprintf('Non-follower load is considered');
    end
    fprintf('\n');
    if ~isscalar(xib) || ~isscalar(etab)
        error('Load extensions xib and etab must be scalar');
    end
    fprintf('Xi parametric coordinate = %.2d\n',xib);
    fprintf('Eta parametric coordinate = %.2d\n',etab);
    fprintf('_______________________________________________________\n\n');

    % start measuring computational time
    tic;
end

%% 0. Read input

% Read the data from the B-Spline patch
p = BSplinePatch.p;
q = BSplinePatch.q;
Xi = BSplinePatch.Xi;
Eta = BSplinePatch.Eta;
CP = BSplinePatch.CP;
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));
if isFollower
    CPd = BSplinePatch.CPd;
end
isNURBS = BSplinePatch.isNURBS;

% Number of Control Points affecting each element
noCPsEl = (p + 1)*(q + 1);

% Number of DOFs per Control Point
noDOFsCP = 3;

% Number of DOFs per element
noDOFsEl = noDOFsCP*noCPsEl;

% Number of DOFs for the patch
if isfield(BSplinePatch,'noDOFs')
    noDOFs = BSplinePatch.noDOFs;
else
    noDOFs = noDOFsCP*nxi*neta;
end

% Initialize auxiliary arrays
RMtx = zeros(noDOFsCP,noDOFsEl);

% Initialize output arrays
F = zeros(noDOFs,1);
if isFollower
    tanMtx = zeros(noDOFs);
else
    tanMtx = 'undefined';
end

%% 1. Find the knot span of the parametric location where the load is applied
xiKnotSpan = findKnotSpan(xib,Xi,nxi);
etaKnotSpan = findKnotSpan(etab,Xi,neta);

%% 2. Get the Element Freedom Table (EFT) corresponding to the element where the load application parametric location belongs
idElmnt = BSplinePatch.knotSpan2ElmntNo(xiKnotSpan,etaKnotSpan);
EFT = BSplinePatch.EFT(:,idElmnt);

%% 3. Compute the IGA basis functions and their derivatives at the parametric location where the load is applied
if strcmp(direction,'theta1') || strcmp(direction,'theta2') || strcmp(direction,'normal')
    noDrv = 1;
else
    noDrv = 0;
end
dR = computeIGABasisFunctionsAndDerivativesForSurface...
    (xiKnotSpan,p,xib,Xi,etaKnotSpan,q,etab,Eta,CP,isNURBS,noDrv);

%% 4. Compute the base vectors
noDrv = 0;
if strcmp(direction,'theta1') || strcmp(direction,'theta2') || strcmp(direction,'normal')
    if isFollower
        [g1,g2] = computeBaseVectorsAndDerivativesForBSplineSurface...
            (xiKnotSpan,p,etaKnotSpan,q,CPd,noDrv,dR);
    else
        [G1,G2] = computeBaseVectorsAndDerivativesForBSplineSurface...
            (xiKnotSpan,p,etaKnotSpan,q,CP,noDrv,dR);
    end
end

%% 5. Compute the normal to the surface vector
if strcmp(direction,'theta1') || strcmp(direction,'theta2') || strcmp(direction,'normal')
    if isFollower
        nTilde = cross(g1,g2);
        n = nTilde/norm(nTilde) ;
    else
        NTilde = cross(G1,G2);
        N = NTilde/norm(NTilde);
    end
end

%% 6. Compute the traction vector

% Compute the amplitude of the force if it is spatially dependent
if ~isnumeric(FAmp)
    X = computeCartesianCoordinatesOfAPointOnBSplineSurface...
        (xiKnotSpan,p,xib,Xi,etaKnotSpan,q,etab,Eta,CP,dR(:,1));
    if  ischar(FAmp) && ~isa(FAmp,'function_handle')
        FAmp = str2func(FAmp);
    end
    FAmplitude = FAmp(X(1,1),X(2,1),X(3,1),t);
else
    FAmplitude = FAmp;
end

% Decide upon the direction
if strcmp(direction,'x')
    dirVct = [1 0 0]';
elseif strcmp(direction,'y')
    dirVct = [0 1 0]';
elseif strcmp(direction,'z')
    dirVct = [0 0 1]';
elseif strcmp(direction,'theta1')
    if isFollower
        dirVct = g1/norm(g1);
    else
        dirVct = G1/norm(G1);
    end
elseif strcmp(direction,'theta2')
    if isFollower
        dirVct = g2/norm(g2);
    else
        dirVct = G2/norm(G2);
    end
elseif strcmp(direction,'normal')
    if isFollower
        dirVct = n;
    else
        dirVct = N;
    end
end

% Compute the tractionvector
tractionVct = FAmplitude*dirVct;

%% 7. Compute the matric containing the basis functions
for iCPs = 1:noCPsEl
    RMtx(1,3*iCPs - 2) = dR(iCPs,1);
    RMtx(2,3*iCPs - 1) = dR(iCPs,1);
    RMtx(3,3*iCPs) = dR(iCPs,1);
end

%% 8. Compute the consistent load vector corresponding to the point load application
F(EFT) = RMtx'*tractionVct;

%% 9. Update by the existing load vector
if isvector(FOutdated)  
    F = F + FOutdated;   
end

%% 10. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    if isvector(FOutdated)
        fprintf('Load vector update took %.2d seconds \n\n',computationalTime);
        fprintf('________________Load Vector Update Ended_______________\n');
        fprintf('#######################################################\n\n\n');
    else
        fprintf('Load vector computation took %.2d seconds \n\n',computationalTime);
        fprintf('____________Load Vector Computation Ended______________\n');
        fprintf('#######################################################\n\n\n');
    end
end

end
