function [F, tanMtxLoad] = ...
    computeLoadVctLinePressureVectorForIGATimoshenkoBeam2D ...
    (Fold, BSplinePatch, xib, etab, loadFnc, loadDir, isFollower, t, ...
    int, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the consistent load vector for the application of line pressure
% load on the isogeometric beam elements of the Timoshenko type.
%
%        Input :
%         Fold : The outdated force vector
% BSplinePatch : Structure containing all computational information for the
%                B-Spline patch,
%                   .p : Polynomial degree of the curve
%                  .Xi : Knot vector of the curve
%                  .CP : Control point coordinates and weights of the curve
%             .isNURBS : Flag on whether the underlying basis is a NURBS or
%                        a  B-Spline
%          xib : Load extension (e.g. [0.5 0.8])
%         etab : Dummy variable for this function
%      loadFnc : Load magnitude or function handle to its computation
%      loadDir : Load direction, 1=tangent 2=normal direction
%   isFollower : Dummy variable for this function
%            t : Time instance
%          int : Gaussian integration parameters
%       outMsg : Whether or not to output message on refinement progress
%                'outputEnabled' : enables output information
%
%       Output :
%            F : Load vector
%   tanMtxLoad : Dummy output for this function
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all elements
% 
%    1i. Loop over all the quadrature points
%
%        1i.1. Transformation from unit interval to the knot span
%
%        1i.2. Compute the NURBS basis functions and their derivatives
%
%        1i.3. Compute the base vectors on the curve
%
%        1i.4. Compute the length of the normal vector and the normal to the 2D curve vector
%
%        1i.5. Compute and assemble the components of the global load vector
%        
% 2. If the given old load vector is not null, sum them up with the newly computed one
%
% 3. Appendix
%
%% Function main body
if strcmp(outMsg, 'outputEnabled')
    fprintf('_____________________________________________________\n');
    fprintf('#####################################################\n');
    fprintf('Computing the consistent load vector corresponding to \n');
    if isnumeric(loadFnc)
        fprintf('constant pressure load on the isogeometric Timoshenko ');
    else
        fprintf('varying pressure load on the isogeometric Timoshenko');
    end
    fprintf('\nbeam reference geometry \n\n');
    fprintf('Pressure application extension equals to [%d,%d]\n', xib(1), xib(2));
    if isnumeric(loadFnc)
        fprintf('Pressure magnitude p = %d \n', loadFnc)
    end
    fprintf('_____________________________________________________\n\n');
    tic;
end

%% 0. Read input

% Get the computational information from the patch
p = BSplinePatch.p;
Xi = BSplinePatch.Xi;
CP = BSplinePatch.CP;
isNURBS = BSplinePatch.isNURBS;

% Initialize the load vector
F = zeros(size(CP, 1)*3, 1);

% Initialize dummy output
tanMtxLoad = 'undefined';

% Number of Control Points
nxi = length(CP(:, 1));

% choose number of gauss points
if strcmp(int.type, 'default')
    ngauss = ceil (0.5*(p + 1));
elseif strcmp(int.type, 'user')
    ngauss = int.ngaussLoad;
end

% Get start end end point for th integration
j1 = findKnotSpan(xib(1), Xi, nxi);
j2 = findKnotSpan(xib(2), Xi, nxi);

% Get the Gauss Points and their coordinates
[GP, GW] = getGaussPointsAndWeightsOverUnitDomain(ngauss);

%% 1. Loop over all elements
for i = j1:j2
    % check if we are in a non-zero element
    if Xi(i + 1) - Xi(i) ~= 0
        %% 1i. Loop over all the quadrature points
        for j = 1:length(GP)     
            %% 1i.1. Transformation from unit interval to the knot span
            xi = ((1 - GP(j))*Xi(i) + (1 + GP(j))*Xi(i + 1))/2;
            xiSpan = findKnotSpan(xi, Xi, nxi);  

            %% 1i.2. Compute the NURBS basis functions and their derivatives
            numDrv = 2;
            dR = computeIGABasisFunctionsAndDerivativesForCurve ...
                (xiSpan, p, xi, Xi, CP, isNURBS, numDrv);

            %% 1i.3. Compute the base vectors on the curve
            [A1, A2] = computeBaseVectorAndNormalToNURBSCurve2D ...
                (xiSpan, p, CP, dR);

            %% 1i.4. Compute the length of the normal vector and the normal to the 2D curve vector
            L = norm(A1);
            n = A2/norm(A2);
            t = A1/L;

            %% 1i.5. Compute and assemble the components of the global load vector
            for v = 1:p + 1
                index1 = 1 + (xiSpan - p - 1)*2 + (v - 1)*2;
                index2 = 2 + (xiSpan - p - 1)*2 + (v - 1)*2;
                indexTim1 = index1 + (xiSpan - p - 1) + v - 1;
                indexTim2 = index2 + (xiSpan - p - 1) + v - 1;
                if loadDir == 1
                    F(index1, 1) = F(index1,1) + ...
                        loadFnc*dR(v, 1)*t(1, 1)*(Xi(i + 1) - Xi(i))/2.0*L*GW(j); 
                    F(index2, 1) = F(index2, 1) + ...
                        loadFnc*dR(v, 1)*t(2, 1)*(Xi(i + 1) - Xi(i))/2.0*L*GW(j); 
                elseif loadDir == 2
                    F(indexTim1, 1) = F(indexTim1, 1) + ...
                        loadFnc*dR(v, 1)*n(1, 1)*(Xi(i + 1) - Xi(i))/2.0*L*GW(j); 
                    F(indexTim2, 1) = F(indexTim2, 1) + ...
                        loadFnc*dR(v, 1)*n(2, 1)*(Xi(i + 1) - Xi(i))/2.0*L*GW(j); 
                else
                    error('Direction %d is not implemented for the loading of a 2D isogeometric Timoshenko beam', loadDir);
                end
            end
        end
    end
end

%% 2. If the given old load vector is not null, sum them up with the newly computed one
if isvector(Fold)  
    F = F + Fold;   
end

%% 3. Appendix
if strcmp(outMsg, 'outputEnabled')
    computationalTime = toc;
    fprintf('Computing the load vector took %d seconds \n\n', computationalTime);
    fprintf('_____________Computing Load Vector Ended_____________\n');
    fprintf('#####################################################\n\n\n');
end

end
