function F = computeLoadVctPointPressureVectorForIGABernoulliBeam2D ...
    (Fold, xib, p, Xi, CP, isNURBS, loadFnc, loadDir, t, int, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the consistent load vector for the application of a point load on 
% the isogeometric beam elements of the Benroulli type.
%
%    Input :
%        p : The polynomial degree of the curve
%      xib : Load extension (e.g. 0.5)
%       Xi : The knot vector of the NURBS curve
%       CP : The Control Points of the NURBS curve
%  isNURBS : Flag on whether the geometrical basis is a NURBS or a B-Spline
%  loadFnc : Pressure magnitude
%  loadDir : Load direction, 1=tangent 2=normal direction
%        t : Time instance
%      int : Gaussian integration parameters
%   outMsg : Whether or not to output message on refinement progress
%            'outputEnabled' : enables output information
%
% Output :
%      F : Load vector 
%
% Function layout :
%
% 0. Read input
%
% 1. Find the knot span
%
% 2. Compute the NURBS basis functions and their derivatives at the load application point
%
% 3. Compute the base vectors on the curve at the load application point
%
% 4. Compute the length of the normal vector and the normal to the 2D curve vector
%
% 5. Compute and assemble the components of the global load vector
%
% 6. If the given old load vector is not null, sum them up with the newly computed one
%
% 7. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('_____________________________________________________\n');
    fprintf('#####################################################\n');
    fprintf('Computing the consistent load vector corresponding to \n');
    if isnumeric(loadFnc)
        fprintf('constant point load on the isogeometric Bernoulli ');
    else
        fprintf('varying point load on the isogeometric Bernoulli ');
    end
    fprintf('\nbeam reference geometry \n\n');
    fprintf('Point load application at xi = %d,%d]\n',xib);
    if isnumeric(loadFnc)
        fprintf('Pressure magnitude p = %d \n',loadFnc)
    end
	fprintf('_____________________________________________________\n\n');
    
    % start measuring computational time
    tic;
end

%% 0. Read input

% Initialize the load vector
F = zeros(size(CP,1)*2,1);

% Number of Control Points
nxi = length(CP(:,1));

%% 1. Find the knot span
xiSpan = findKnotSpan(xib,Xi,nxi);

%% 2. Compute the NURBS basis functions and their derivatives at the load application point
nDeriv = 2;
dR = computeIGABasisFunctionsAndDerivativesForCurve...
    (xiSpan,p,xib,Xi,CP,isNURBS,nDeriv);

%% 3. Compute the base vectors on the curve at the load application point
[A1,A2] = computeBaseVectorAndNormalToNURBSCurve2D...
    (xiSpan,p,CP,dR);

%% 4. Compute the length of the normal vector and the normal to the 2D curve vector
L = norm(A1);
n = A2/norm(A2);
t = A1/L;

%% 5. Compute and assemble the components of the global load vector
for v = 1:p+1
    index1 = 1 + (xiSpan-p-1)*2 + (v-1)*2;
    index2 = 2 + (xiSpan-p-1)*2 + (v-1)*2;
    if loadDir == 1
        F(index1,1) = F(index1,1) + loadFnc*dR(v,1)*t(1,1); 
        F(index2,1) = F(index2,1) + loadFnc*dR(v,1)*t(2,1); 
    elseif loadDir == 2
        F(index1,1) = F(index1,1) + loadFnc*dR(v,1)*n(1,1); 
        F(index2,1) = F(index2,1) + loadFnc*dR(v,1)*n(2,1);
    end
end

%% 6. If the given old load vector is not null, sum them up with the newly computed one
if isvector(Fold)  
    F = F + Fold;   
end

%% 7. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    fprintf('Computing the load vector took %d seconds \n\n',computationalTime);
    fprintf('_____________Computing Load Vector Ended_____________\n');
    fprintf('#####################################################\n\n\n');
end

end
