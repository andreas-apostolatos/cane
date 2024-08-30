function [F, tanMtx] = computeLoadVctLineIGAPlateInMembraneAction ...
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
% Returns the consistent nodal forces to a line load fload [N/m] for the 
% isogeometric plate in membrane action problem. The direction of FAmp can 
% be in x, y, parallel or perpendicular to the edge.
%
%          Input :
%      FOutdated : The existing force vector
%   BSplinePatch : Array containing information on the B-Spline patch
%                      .p,.q : Polynomial orders in xi- and -eta direction
%                   .Xi,.Eta : Polynomial orders in xi-eta direction
%                        .CP : Set of Control Point coordinates and weights
%                   .isNURBS : Flag on whether the basis is a NURBS or a
%                              B-Spline
%                    .noDOFs : Total Number of DOFs for the B-Spline patch
%                              including possible Lagrange Multipliers
%       xib,etab : load extension (e.g. xib = [0 .7], etab = [.5 1])
%           FAmp : The amplitude of the constant line load or handle to
%                  load function [N/m] for varying loads
%      direction : Direction of the applied force
%                           'x' : along x -> FAmp*e1 
%                           'y' : along y -> FAmp*e2
%                           'z' : along z -> FAmp*e3
%                      'theta1' : along theta1 -> FAmp*G1/norm(G1)
%                      'theta2' : along theta2 -> FAmp*G2/norm(G2)
%                      'normal' : normal to the surface -> FAmp*n
%     isFollower : Flag on whether the applied load is a follower load
%              t : The time instance
%            int : On the quadrature,
%                           .type : 'user' or 'default'
%                   .xiNGPForLoad : Number of Gauss Points along xi
%                  .etaNGPForLoad : Number of Gauss Points along eta
%         outMsg : Whether or not to output message on the load vector
%                  computation,
%                   'outputEnabled' : enables output information
%
%         Output :
%              F : The updated load vector
%         tanMtx : Dummy output for this function
%
% Function layout :
%
% 0. Read input
%
% 1. Get the parametrization of the line load for the integration
%
% 2. Get a quadrature rule
%
% 3. Loop over all elements of the load application
% ->
%    3i. Compute the determinant of the Jacobian corresponding to the transformation from the parameter to the integration space
%
%   3ii. Get the Element Freedom Table for the given element
%
%  3iii. Loop over all the Gauss Points
%  ->
%        3iii.1. Compute the parametric coordinate of the Gauss Point
%
%        3iii.2. Compute the IGA basis functions at the Gauss Point
%
%        3iii.3. Compute the base vectors and their derivatives
%
%        3iii.4. Compute the determinant of the Jacobian of the transformation from the Cartesian to the paramter space
%
%        3iii.5. Compute the traction vector
%
%        3iii.6. Compute the matrices containing the basis functions and their parametric derivatives
%
%        3iii.7. Compute the elementary area content at the Gauss Point
%
%        3iii.8. Compute and assemble the element load vector and the tangent matrix contribution at the Gauss Point
%  <-
% <-
%
% 4. Update the existing load vector
%
% 5. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('_______________________________________________________\n');
    fprintf('#######################################################\n');
    if isvector(FOutdated)
        fprintf('Update of the load vector corresponding to line\n');
    else
        fprintf('Computation of the load vector corresponding to line\n');
    end
    fprintf('boundary load for the isogeometric Kirchhoff-Love shell\n');
    fprintf('problem has been initiated\n\n');
    if isnumeric(FAmp)
        fprintf('Constant boundary load is assumed with amplitude = %.2d\n', FAmp);
    else
        fprintf('Varying boundary load is assumed\n');
    end
    if strcmp(direction, 'x') || strcmp(direction, 'y') ...
            || strcmp(direction, 'theta1') || strcmp(direction, 'theta2')
        fprintf('Load direction is chosen as %s', direction);
    else
        error('Select a valid direction for the load');
    end
    fprintf('\n');
    if isFollower
        fprintf('Follower load for this type analysis has not been implemented');
    else
        fprintf('Non-follower load is considered');
    end
    fprintf('\n');
    if isscalar(xib) && isscalar(etab)
        error('Both load extensions xib and etab cannot have simultaneously scalar values for a line load');
    end
    if ~isscalar(xib) && ~isscalar(etab)
        error('Both load extensions xib and etab cannot have simultaneously vactor values for a line load');
    end
    if isscalar(xib)
        fprintf('Load extension in xi parametric direction = [%.2d,%.2d]\n', xib, xib);
    else
        fprintf('Load extension in xi parametric direction = [%.2d,%.2d]\n', xib(1), xib(2));
    end
    if isscalar(etab)
        fprintf('Load extension in eta parametric direction = [%.2d,%.2d]\n', etab, etab);
    else
        fprintf('Load extension in eta parametric direction = [%.2d,%.2d]\n', etab(1), etab(2));
    end
    fprintf('_______________________________________________________\n\n');
    tic;
end

%% 0. Read input

% Read the data from the B-Spline patch
p = BSplinePatch.p;
q = BSplinePatch.q;
Xi = BSplinePatch.Xi;
Eta = BSplinePatch.Eta;
CP = BSplinePatch.CP;
nxi = length(CP(:, 1, 1));
neta = length(CP(1, :, 1));
isNURBS = BSplinePatch.isNURBS;

% Number of Control Points affecting each element
numCPsEl = (p + 1)*(q + 1);

% Number of DOFs per Control Point
numDOFsCP = 2;

% Number of DOFs per element
numDOFsEl = numDOFsCP*numCPsEl;

% Number of derivatives for the B-Spline basis functions
numDrvsBasis = 1;
numDrvsBaseVct = 0;

% Number of DOFs for the patch
if isfield(BSplinePatch, 'noDOFs')
    numDOFs = BSplinePatch.noDOFs;
else
    numDOFs = numDOFsCP*nxi*neta;
end

% Initialize auxiliary arrays
RMtx = zeros(numDOFsCP, numDOFsEl);

% Initialize output arrays
F = zeros(numDOFs, 1);
tanMtx = 'undefined';

%% 1. Get the parametrization of the line load for the integration
if isscalar(xib)
    % Set the flag on whether the integration takes place along xi to false
    isOnXi = false;
    
    % Get the fixed parameter
    xi = xib;
    
    % Find the knot span index of the fixed parameter
    xiKnotSpan = findKnotSpan(xi, Xi, nxi);
    
    % Get the knot vector over which the integration takes place
    thetaKnotVct = Eta; 
    
    % Set the starting and ending parametric coordinates for the 
    % integration
    thetaStart = etab(1);
    thetaStartKnotSpan = findKnotSpan(thetaStart, Eta, neta);
    thetaEnd = etab(2);
    thetaEndKnotSpan = findKnotSpan(thetaEnd, Eta, neta);
else
    % Set the flag on whether the integration takes place along xi to true
    isOnXi = true;
    
    % Get the fixed parameter
    eta = etab;
    
    % Find the knot span index of the fixed parameter
    etaKnotSpan = findKnotSpan(eta, Eta, neta);
    
    % Get the knot vector over which the integration takes place
    thetaKnotVct = Xi; 
    
    % Set the starting and ending parametric coordinates for the 
    % integration
    thetaStart = xib(1);
    thetaStartKnotSpan = findKnotSpan(thetaStart, Xi, nxi);
    thetaEnd = xib(2);
    thetaEndKnotSpan = findKnotSpan(thetaEnd, Xi, nxi);
end

%% 2. Get a quadrature rule
if isstruct(propInt)
    if ~isfield(propInt, 'type')
        error('propInt must define a type')
    end
    if strcmp(propInt.type, 'default')
        if isOnXi
            numGPs = ceil(2*p - 1);
        else
            numGPs = ceil(2*q - 1);
        end
    elseif strcmp(propInt.type, 'user')
        if isOnXi
            numGPs = propInt.xiNGPForLoad;
        else
            numGPs = propInt.etaNGPForLoad;
        end
    else
        error('Select a valid integration type int.type')
    end
else
    error('propInt must be a structure');
end
[GP, GW] = getGaussPointsAndWeightsOverUnitDomain(numGPs);

%% 3. Loop over all elements of the load application
for iKnotSpan = thetaStartKnotSpan:thetaEndKnotSpan
    if thetaKnotVct(iKnotSpan + 1) ~= thetaKnotVct(iKnotSpan)
        %% 3i. Compute the determinant of the Jacobian corresponding to the transformation from the parameter to the integration space
        detJParam2Integr = (thetaKnotVct(iKnotSpan + 1) - thetaKnotVct(iKnotSpan))/2;
        
        %% 3ii. Get the Element Freedom Table for the given element
        if isOnXi
            xiKnotSpan = iKnotSpan;
        else
            etaKnotSpan = iKnotSpan;
        end
        noElmnt = BSplinePatch.knotSpan2ElmntNo(xiKnotSpan, etaKnotSpan);
        EFT = BSplinePatch.EFT(:, noElmnt);
        
        %% 3iii. Loop over all the Gauss Points
        for iGP = 1:numGPs
            %% 3iii.1. Compute the parametric coordinate of the Gauss Point
            theta = (thetaKnotVct(iKnotSpan + 1) + thetaKnotVct(iKnotSpan) + ...
                GP(iGP)*(thetaKnotVct(iKnotSpan + 1) - thetaKnotVct(iKnotSpan)))/2;
            
            %% 3iii.2. Compute the IGA basis functions at the Gauss Point
            if isOnXi
                xi = theta;
                xiKnotSpan = findKnotSpan(xi, Xi, nxi);
            else
                eta = theta;
                etaKnotSpan = findKnotSpan(eta, Eta, neta);
            end
            dR = computeIGABasisFunctionsAndDerivativesForSurface ...
                (xiKnotSpan, p, xi, Xi, etaKnotSpan, q, eta, Eta, CP, ...
                isNURBS, numDrvsBasis);
            
            %% 3iii.3. Compute the base vectors and their derivatives
            [A1, A2] = computeBaseVectorsAndDerivativesForBSplineSurface ...
                (xiKnotSpan, p, etaKnotSpan, q, CP, numDrvsBaseVct, dR);
            
            %% 3iii.4. Compute the determinant of the Jacobian of the transformation from the Cartesian to the paramter space
            if isOnXi
                detJCartesian2Param = norm(A1);
            else
                detJCartesian2Param = norm(A2);
            end
            
            %% 3iii.5. Compute the traction vector
                    
            % Compute the amplitude of the force if it is spatially 
            % dependent
            if ~isnumeric(FAmp)
                X = computeCartesianCoordinatesOfAPointOnBSplineSurface ...
                    (xiKnotSpan, p, xi, Xi, etaKnotSpan, q, eta, Eta, ...
                    CP, dR(:, 1));
                if  ischar(FAmp) && ~isa(FAmp, 'function_handle')
                    FAmp = str2func(FAmp);
                end
                FAmplitude = FAmp(X(1, 1), X(2, 1), X(3, 1), t);
            else
                FAmplitude = FAmp;
            end

            % Decide upon the direction
            if strcmp(direction, 'x')
                dirVct = [1; 0];
            elseif strcmp(direction, 'y')
                dirVct = [0; 1];
            elseif strcmp(direction,'theta1')
                dirVct = A1(:, 1)/norm(A1(:, 1));
            elseif strcmp(direction, 'theta2')
                dirVct = A2(:, 1)/norm(A2(:, 1));
            end

            % Compute the tractionvector
            tractionVct = FAmplitude*dirVct;
            
            %% 3iii.6. Compute the matrices containing the basis functions and their parametric derivatives
            for iCPs = 1:numCPsEl
                RMtx(1, 2*iCPs - 1) = dR(iCPs,1);
                RMtx(2, 2*iCPs) = dR(iCPs,1);
            end
            
            %% 3iii.7. Compute the elementary area content at the Gauss Point
            elementaryAreaGP = detJCartesian2Param*detJParam2Integr*GW(iGP);
            
            %% 3iii.8. Compute and assemble the element load vector and the tangent matrix contribution at the Gauss Point
            F(EFT) = F(EFT) + RMtx'*tractionVct*elementaryAreaGP;
        end
    end
end

%% 4. Update the existing load vector
if isvector(FOutdated)  
    F = F + FOutdated;
end

%% 5. Appendix
if strcmp(outMsg, 'outputEnabled')
    computationalTime = toc;
    if isvector(FOutdated)
        fprintf('Load vector update took %.2d seconds \n\n', computationalTime);
        fprintf('______________Load Vector Update Ended______________\n');
        fprintf('####################################################\n\n\n');
    else
        fprintf('Load vector computation took %.2d seconds \n\n', computationalTime);
        fprintf('____________Load Vector Computation Ended___________\n');
        fprintf('####################################################\n\n\n');
    end
end

end