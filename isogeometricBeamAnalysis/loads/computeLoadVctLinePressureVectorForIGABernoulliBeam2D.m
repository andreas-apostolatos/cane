function F = computeLoadVctLinePressureVectorForIGABernoulliBeam2D ...
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
% Returns the consistent load vector for the application of line pressure
% load on the isogeometric beam elements of the Benroulli type.
%
%    Input :
%        p : The polynomial degree of the curve
%      xib : Load extension (e.g. [0.5 0.8])
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
if strcmp(outMsg,'outputEnabled')
    fprintf('_____________________________________________________\n');
    fprintf('#####################################################\n');
    fprintf('Computing the consistent load vector corresponding to \n');
    if isnumeric(loadFnc)
        fprintf('constant pressure load on the isogeometric Bernoulli ');
    else
        fprintf('varying pressure load on the isogeometric Bernoulli ');
    end
    fprintf('\nbeam reference geometry \n\n');
    fprintf('Pressure application extension equals to [%d,%d]\n',xib(1),xib(2));
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

% choose number of gauss points
if strcmp(int.type,'default')
    noGP = ceil ( 0.5*(p + 1));
elseif strcmp(int.type,'user')
    noGP = int.noGPLoad;
end

% Get start end end point for th integration
j1 = findKnotSpan(xib(1),Xi,nxi);
j2 = findKnotSpan(xib(2),Xi,nxi);

% Get the Gauss Points and their coordinates
[GP,GW] = getGaussPointsAndWeightsOverUnitDomain(noGP);

%% 1. Loop over all elements
for i=j1:j2
    % check if we are in a non-zero element
    if Xi(i+1)-Xi(i)~=0  
        %% 1i. Loop over all the quadrature points
        for j = 1 : length(GP)     
            %% 1i.1. Transformation from unit interval to the knot span
            xi = ((1-GP(j))*Xi(i)+(1+GP(j))*Xi(i+1))/2;                                  
            xiSpan = findKnotSpan(xi,Xi,nxi);  

            %% 1i.2. Compute the NURBS basis functions and their derivatives
            nDeriv = 2;
            dR = computeIGABasisFunctionsAndDerivativesForCurve...
                (xiSpan,p,xi,Xi,CP,isNURBS,nDeriv);

            %% 1i.3. Compute the base vectors on the curve
            [A1,A2] = computeBaseVectorAndNormalToNURBSCurve2D...
                (xiSpan,p,CP,dR);

            %% 1i.4. Compute the length of the normal vector and the normal to the 2D curve vector
            L = norm(A1);
            n = A2/norm(A2);
            t = A1/L;

            %% 1i.5. Compute and assemble the components of the global load vector
            for v=1:p+1
                index1 = 1 + (xiSpan-p-1)*2 + (v-1)*2;
                index2 = 2 + (xiSpan-p-1)*2 + (v-1)*2;
                if loadDir == 1
                    F(index1,1) = F(index1,1) + loadFnc*dR(v,1)*t(1,1)*(Xi(i+1)-Xi(i))/2.0*L*GW(j); 
                    F(index2,1) = F(index2,1) + loadFnc*dR(v,1)*t(2,1)*(Xi(i+1)-Xi(i))/2.0*L*GW(j); 
                elseif loadDir == 2
                    F(index1,1) = F(index1,1) + loadFnc*dR(v,1)*n(1,1)*(Xi(i+1)-Xi(i))/2.0*L*GW(j); 
                    F(index2,1) = F(index2,1) + loadFnc*dR(v,1)*n(2,1)*(Xi(i+1)-Xi(i))/2.0*L*GW(j);
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
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    fprintf('Computing the load vector took %d seconds \n\n',computationalTime);
    fprintf('_____________Computing Load Vector Ended_____________\n');
    fprintf('#####################################################\n\n\n');
end

end
