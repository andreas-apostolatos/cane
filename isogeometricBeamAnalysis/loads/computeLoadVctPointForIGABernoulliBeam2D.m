function Fl = computeLoadVctPointForIGABernoulliBeam2D ...
    (Fold, xib, p, Xi, CP, isNURBS, loadAmp, loadDir, t, int, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the consistent load vector corresponding to the application of a
% point load corresponding to the isogeometric beam of Bernoulli type.
%
%    Input :
%     Fold : The outdated load vector
%        p : The polynomial degree of the curve
%      xib : Parametric location on where the load is applied
%       Xi : The knot vector of the NURBS curve
%       CP : The Control Points of the NURBS curve
%  isNURBS : Flag on whether the geometrical basis is a NURBS or a B-Spline
%  loadAmp : Load magnitude
%  loadDir : Load direction, 1=X 2=Y direction
%        t : Time instance
%      int : Gaussian integration parameters
%   outMsg : Whether or not to output message on refinement progress
%            'outputEnabled' : enables output information
%
% Output :
%     Fl : the updated Load vector
%
% Function layout :
%
% 0. Read input
%
% 1. Compute the consistent load vector corresponding to a point load
%
% 2. Re-assemble the load vector with respect to the global numbering
%
% 3. If the given old load vector is not null, sum them up with the newly computed one
%
% 4. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('_____________________________________________________\n');
    fprintf('#####################################################\n');
    fprintf('Computing the consistent load vector corresponding to \n');
    fprintf('constant pressure load on the isogeometric Bernoulli ');
    fprintf('\nbeam reference geometry \n\n');
    if ~isscalar(xib)
       error('Wrong load extension in load application');
    end
    fprintf('Load is applied at xi = %d\n',xib);
    fprintf('Load direction is dir = %d\n',loadDir);
	fprintf('_____________________________________________________\n\n');
    
    % start measuring computational time
    tic;
end

%% 0. Read input

% Number of Control Points
nxi = length(CP(:,1));

% Number of Control Points at the element level
noCPsEl = (p + 1);

% Number of DOFs at the global level
noDOFs = 2*nxi;

% Initialize the element load vector
FEl = zeros(noCPsEl,1);

% Initialize the load vector
F = zeros(nxi,2);

% Initialize output array
Fl = zeros(noDOFs,1);

%% 1. Compute the consistent load vector corresponding to a point load
xiSpan = findKnotSpan(xib,Xi,nxi);
R = computeIGABasisFunctionsAndDerivativesForCurve...
    (xiSpan,p,xib,Xi,CP,isNURBS,0);
for counterBasis = 1:p+1
    FEl(counterBasis,1) = loadAmp*R(counterBasis,1); 
end
F(xiSpan - p:xiSpan,loadDir) = FEl;

%% 2. Re-assemble the load vector with respect to the global numbering
counter = 1;
for i = 1:length(F(:,1))
    % Assemble the x-coordinates of the load vector
    Fl(counter,1,1) = F(i,1);

    % Assemble the y-coordinates of the load vector
    Fl(counter + 1,1) = F(i,2);

    % Update counter
    counter = counter + 2;
end

%% 3. If the given old load vector is not null, sum them up with the newly computed one
if isvector(Fold)  
    Fl = Fl + Fold;   
end

%% 4. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    fprintf('Computing the load vector took %d seconds \n\n',computationalTime);
    fprintf('_____________Computing Load Vector Ended_____________\n');
    fprintf('#####################################################\n\n\n');
end

end
