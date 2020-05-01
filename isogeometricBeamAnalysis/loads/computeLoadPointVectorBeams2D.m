function F = computeLoadPointVectorBeams2D...
    (Fold, p, Xi, CP, isNURBS, loadAmplitude, xib, analysis, dir)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% returns the consistent nodal forces to a point load loadAmplitude for the 
% isogeometric 2D beams of type Bernoulli and Timoshenko.
%
%         Input :
%          Fold : existing force vector
%           xib : load position
%             p : polynomial degree in u,v-direction
%            Xi : Knot vectors in u,v-direction
%            CP : Set of Control Point coordinates and weights
%       isNURBS : Flag on whether the basis is a NURBS or a B-Spline
% loadAmplitude : point load magnitude
%      analysis : .type : Bernoulli or Timoshenko
%           dir : direction of f  1=x, 2=y
%
%        Output :
%             F : updated force vector
%
% 0. Read input
%
% 1. Evaluate the NURBS basis functions at the point of the load application
%
% 2. Compute the load component at the application point
%
% 3. If the given old load vector is not null, sum them up with the newly computed one
%
%% Function main body

%% 0. Read input
nu = length(CP(:,1));

i = findKnotSpan(xib,Xi,nu);

% Initialize the load vector
if strcmp(analysis.type,'Bernoulli')
    F = zeros(size(CP,1)*2,1);
elseif strcmp(analysis.type,'Timoshenko')
    F = zeros(size(CP,1)*3,1);
end

%% 1. Evaluate the NURBS basis functions at the point of the load application
R = computeIGABasisFunctionsAndDerivativesForCurve(i,p,xib,Xi,CP,isNURBS,0);

%% 2. Compute the load component at the application point

% Initialize counter
k = 1;

% Loop over all the contributions at the current knot span
    
for b=1:p+1
    % Find the correct index for the current knot span
    index = dir + (i-p-1)*2 + (b-1)*2;
    
    % Compute the entries of the load vector according to the respective
    % theory
    if  strcmp(analysis.type,'Bernoulli')
        F(index,1) = R(k)*loadAmplitude; 
    elseif strcmp(analysis.type,'Timoshenko')
        index_tim = index + (i-p-1) + b - 1;
        F(index_tim,1) = R(k)*loadAmplitude; 
    end

    % Update counter
    k = k + 1;

end
    
%% 3. If the given old load vector is not null, sum them up with the newly computed one
if isvector(Fold)  
    F = F + Fold;   
end

end
