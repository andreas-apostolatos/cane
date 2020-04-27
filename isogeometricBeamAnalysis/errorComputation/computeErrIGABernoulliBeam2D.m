function [relErrL2, minElSize] = computeErrIGABernoulliBeam2D...
    (p, Xi, CP, isNURBS, dHat,computeExactSolution, problemSettings, ...
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
% Returns the relative error in the L2-norm for the case of an isogeometric
% Bernoulli beam element given the rule for the computation of the exact
% solution.
%
%                Input :
%                    p : The polynomial order of the B-Spline patch
%                   Xi : The knot vector of the B-Spline patch
%                   CP : The set of Control Point coordinates and weights 
%                        for the B-Spline patch
%              isNURBS : Flag on whether the underlying patch is a NURBS or 
%                        a B-Spline
%                 dHat : The discrete solution vector
% computeExactSolution : Handle for a function which returns the exact
%                        solution for a given Euler-Bernoulli beam problem
%      problemSettings : Parameters which are related to the benchmark
%                        problem, such as dimensions, load amplitude etc.
%              propInt : On the numerical quadrature
%               outMsg : On outputting information
%
%               Output :
%                errL2 : The error in the L2-norm
%            minElSize : The minimum element size in the isogeometric mesh
%
% Function layout :
%
% 0. Read input
%
% 1. Get Gauss Points and Gauss coordinates
%
% 2. Loop over all the elements on the mesh
% ->
%    2i. Initialize the element area size
%
%   2ii. Compute the determinant of the Jacobian to the transformation from the parameter to the integration space
%
%  2iii. Loop over all the Gauss Points
%  ->
%        2iii.1. Linear transformation from the quadrature domain to the knot span
%
%        2iii.2. Find the knot span
%
%        2iii.3. Compute the IGA basis functions and their derivatives
%
%        2iii.4. Compute the determinant of the transformation from the physical to the parameter space
%
%        2iii.5. Compute the displacement field at the Gauss point and the physical coordinates of the Gauss Point
%
%        2iii.6. Get the analytical displacement field on the Gauss Point
%
%        2iii.7. Compute the relative error on the Gauss point and add the contribution
%  <-
%   2iv. Check if the current minimum element area size is smaller than for the previous element
% <-
% 3. Compute the relative error in the L2-norm
%
% 4. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('________________________________________________________\n');
    fprintf('########################################################\n');
    fprintf('Computation the relative error in the displacement field\n');
    fprintf('for the isogeometric Euler-Bernoulli beam has been \n');
    fprintf('initiated\n');
    fprintf('________________________________________________________\n\n');

    % start measuring computational time
    tic;
end

%% 0. Read input

% Number of Control Points
nxi = length(CP(:,1));

% Initialize minimum element area size
minElementSize = norm(CP(1,1:3) - CP(nxi,1:3));
minElSize = minElementSize;

% Initialize auxiliary variables
errL2 = 0;
exactL2 = 0;

%% 1. Get Gauss Points and Gauss coordinates
if strcmp(propInt.type,'default')
    nGP = ceil((p + 1)/2);
elseif strcmp(propInt.type,'user')
    nGP = propInt.noGP;
end
[GP,GW] = getGaussPointsAndWeightsOverUnitDomain(nGP);

%% 2. Loop over all the elements on the mesh
for i = p+1:length(Xi)-p-1
    % Check if we are in a non zero span
    if Xi(i+1)-Xi(i)~=0
        %% 2i. Initialize the element area size
            
        % Initialize the new one
        minElementSize = 0;
        
        %% 2ii. Compute the determinant of the Jacobian to the transformation from the parameter to the integration space
        detJxiu = (Xi(i+1)-Xi(i))/2.0;
        
        %% 2iii. Loop over all the Gauss Points
        for j = 1:length(GW)
            %% 2iii.1. Linear transformation from the quadrature domain to the knot span
            xi = ((1-GP(j))*Xi(i)+(1+GP(j))*Xi(i+1))/2;
            
            %% 2iii.2. Find the knot span
            xiSpan = findKnotSpan(xi,Xi,nxi);
            
            %% 2iii.3. Compute the IGA basis functions and their derivatives
            dR = computeIGABasisFunctionsAndDerivativesForCurve(xiSpan,p,xi,Xi,CP,isNURBS,1);
            
            %% 2iii.4. Compute the determinant of the transformation from the physical to the parameter space
            [G,~] = computeBaseVectorAndNormalToNURBSCurve2D(xiSpan,p,CP,dR);
            detJxxi = norm(G);
            
            %% 2iii.5. Compute the displacement field at the Gauss point and the physical coordinates of the Gauss Point
            dIGA = zeros(2,1);
            P = zeros(3,1);
            for b = 0:p
                % Compute the displacement components iteratively
                dIGA(1,1) = dIGA(1,1) + dR(b+1,1)*dHat(2*(xiSpan-p+b)-1,1);
                dIGA(2,1) = dIGA(2,1) + dR(b+1,1)*dHat(2*(xiSpan-p+b),1);
                
                % Compute the Cartesian coordinates of the Gauss point
                % iteratively
                for a = 1:3
                    P(a,1) = P(a,1) + dR(b+1,1)*CP(xiSpan-p+b,a);
                end
            end
            
            %% 2iii.6. Get the analytical displacement field on the Gauss Point
            dExact = computeExactSolution(P,problemSettings);
            
            %% 2iii.7. Compute the relative error on the Gauss point and add the contribution
            errL2 = errL2 + norm(dIGA - dExact)^2*detJxiu*detJxxi;
            exactL2 = exactL2 + norm(dExact)^2*detJxiu*detJxxi;
        end
    end
    %% 2iv. Check if the current minimum element area size is smaller than for the previous element
    if minElementSize <= minElSize
        minElSize = minElementSize;
    end
end
%% 3. Compute the relative error in the L2-norm
relErrL2 = sqrt(errL2)/sqrt(exactL2);
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Relative displacement error in the L2-norm = %d\n',relErrL2);
    fprintf('>> Minimum element size = %d \n\n',minElSize);
end

%% 4. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    fprintf('Error computation took %d seconds \n\n',computationalTime);
    fprintf('_________________Error Computation Ended________________\n');
    fprintf('########################################################\n\n\n');
end

end
