function [F,tanMtx] = computeLoadVctLineIGAThinStructure...
    (FOutdated,BSplinePatch,xib,etab,FAmp,direction,isFollower,t,int,outMsg)
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
% corresponding to a boundary (line) load.
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
%       xib,etab : load extension (e.g. xib = [0 .7], etab = 1)
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
%         tanMtx : The tangent matrix contribution in case the load is a
%                  follower load
%
% Function layout :
%
% 0. Read input
%
% 1. Get the parametrization of the line load for the integration
%
% 2. Get a quadrature rule
%
% 3. Loop over all elements at the boundary where the load is applied
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
%        3iii.5. Compute the normal to the surface vector
%
%        3iii.6. Compute the traction vector
%
%        3iii.7. Compute the matrices containing the basis functions and their parametric derivatives
%
%        3iii.8. Compute the elementary area content at the Gauss Point
%
%        3iii.9. Compute and assemble the element load vector and the tangent matrix contribution at the Gauss Point
%  <-
% <-
%
% 4. Update by the existing load vector
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
        fprintf('Constant boundary load is assumed with amplitude = %.2d\n',FAmp);
    else
        fprintf('Varying boundary load is assumed\n');
    end
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
    if isscalar(xib) && isscalar(etab)
        error('Both load extensions xib and etab cannot have simultaneously scalar values for a line load');
    end
    if ~isscalar(xib) && ~isscalar(etab)
        error('Both load extensions xib and etab cannot have simultaneously vactor values for a line load');
    end
    if isscalar(xib)
        fprintf('Load extension in xi parametric direction = [%.2d,%.2d]\n',xib,xib);
    else
        fprintf('Load extension in xi parametric direction = [%.2d,%.2d]\n',xib(1),xib(2));
    end
    if isscalar(etab)
        fprintf('Load extension in eta parametric direction = [%.2d,%.2d]\n',etab,etab);
    else
        fprintf('Load extension in eta parametric direction = [%.2d,%.2d]\n',etab(1),etab(2));
    end
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

%% 1. Get the parametrization of the line load for the integration
if isscalar(xib)
    % Set the flag on whether the integration takes place along xi to false
    isOnXi = false;
    
    % Get the fixed parameter
    xi = xib;
    
    % Find the knot span index of the fixed parameter
    xiKnotSpan = findKnotSpan(xi,Xi,nxi);
    
    % Get the knot vector over which the integration takes place
    thetaKnotVct = Eta; 
    
    % Set the starting and ending parametric coordinates for the 
    % integration
    thetaStart = etab(1);
    thetaStartKnotSpan = findKnotSpan(thetaStart,Eta,neta);
    thetaEnd = etab(2);
    thetaEndKnotSpan = findKnotSpan(thetaEnd,Eta,neta);
else
    % Set the flag on whether the integration takes place along xi to true
    isOnXi = true;
    
    % Get the fixed parameter
    eta = etab;
    
    % Find the knot span index of the fixed parameter
    etaKnotSpan = findKnotSpan(eta,Eta,neta);
    
    % Get the knot vector over which the integration takes place
    thetaKnotVct = Xi; 
    
    % Set the starting and ending parametric coordinates for the 
    % integration
    thetaStart = xib(1);
    thetaStartKnotSpan = findKnotSpan(thetaStart,Xi,nxi);
    thetaEnd = xib(2);
    thetaEndKnotSpan = findKnotSpan(thetaEnd,Xi,nxi);
end

%% 2. Get a quadrature rule
if strcmp(int.type,'default')
    if isOnXi
        noGPs = ceil(2*p - 1);
    else
        noGPs = ceil(2*q - 1);
    end
elseif strcmp(int.type,'user')
    if isOnXi
        noGPs = int.xiNGPForLoad;
    else
        noGPs = int.etaNGPForLoad;
    end
else
    error('Select a valid integration type int.type')
end
[GP,GW] = getGaussPointsAndWeightsOverUnitDomain(noGPs);

%% 3. Loop over all elements at the boundary where the load is applied
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
        noElmnt = BSplinePatch.knotSpan2ElmntNo(xiKnotSpan,etaKnotSpan);
        EFT = BSplinePatch.EFT(:,noElmnt);
        
        %% 3iii. Loop over all the Gauss Points
        for iGP = 1:noGPs
            %% 3iii.1. Compute the parametric coordinate of the Gauss Point
            theta = (thetaKnotVct(iKnotSpan + 1) + thetaKnotVct(iKnotSpan) + GP(iGP)*(thetaKnotVct(iKnotSpan + 1) - thetaKnotVct(iKnotSpan)))/2;
            
            %% 3iii.2. Compute the IGA basis functions at the Gauss Point
            if isOnXi
                xi = theta;
                xiKnotSpan = findKnotSpan(xi,Xi,nxi);
            else
                eta = theta;
                etaKnotSpan = findKnotSpan(eta,Eta,neta);
            end
            noDrv = 1;
            dR = computeIGABasisFunctionsAndDerivativesForSurface...
                (xiKnotSpan,p,xi,Xi,etaKnotSpan,q,eta,Eta,CP,isNURBS,noDrv);
            
            %% 3iii.3. Compute the base vectors and their derivatives
            noDrv = 0;
            if isFollower
                [a1,a2] = computeBaseVectorsAndDerivativesForBSplineSurface...
                    (xiKnotSpan,p,etaKnotSpan,q,CPd,noDrv,dR);
            else
                [A1,A2] = computeBaseVectorsAndDerivativesForBSplineSurface...
                    (xiKnotSpan,p,etaKnotSpan,q,CP,noDrv,dR);
            end
            
            %% 3iii.4. Compute the determinant of the Jacobian of the transformation from the Cartesian to the paramter space
            if isOnXi
                if isFollower
                    detJCartesian2Param = norm(a1);
                else
                    detJCartesian2Param = norm(A1);
                end
            else
                if isFollower
                    detJCartesian2Param = norm(a2);
                else
                    detJCartesian2Param = norm(A2);
                end
            end
            
            %% 3iii.5. Compute the normal to the surface vector
            if isFollower
                a3Tilde = cross(a1(:,1),a2(:,1));
                a3 = a3Tilde/norm(a3Tilde) ; 
            else
                A3Tilde = cross(A1(:,1),A2(:,1));
                A3 = A3Tilde/norm(A3Tilde);
            end
            
            %% 3iii.6. Compute the traction vector
                    
            % Compute the amplitude of the force if it is spatially 
            % dependent
            if ~isnumeric(FAmp)
                X = computeCartesianCoordinatesOfAPointOnBSplineSurface...
                    (iXiSpan,p,xi,Xi,iEtaSpan,q,eta,Eta,CP,dR(:,1));
                if  ischar(FAmp) && ~isa(FAmp,'function_handle')
                    FAmp = str2func(FAmp);
                end
                FAmplitude = FAmp(X(1,1),X(2,1),X(3,1),t);
            else
                FAmplitude = FAmp;
            end

            % Decide upon the direction
            if strcmp(direction,'x')
                dirVct = [1; 0; 0];
            elseif strcmp(direction,'y')
                dirVct = [0; 1; 0];
            elseif strcmp(direction,'z')
                dirVct = [0; 0; 1];
            elseif strcmp(direction,'theta1')
                if isFollower
                    dirVct = a1(:,1)/norm(a1(:,1));
                else
                    dirVct = A1(:,1)/norm(A1(:,1));
                end
            elseif strcmp(direction,'theta2')
                if isFollower
                    dirVct = a2(:,1)/norm(a2(:,1));
                else
                    dirVct = A2(:,1)/norm(A2(:,1));
                end
            elseif strcmp(direction,'normal')
                if isFollower
                    dirVct = a3;
                else
                    dirVct = A3;
                end
            end

            % Compute the tractionvector
            tractionVct = FAmplitude*dirVct;
            
            %% 3iii.7. Compute the matrices containing the basis functions and their parametric derivatives
            for iCPs = 1:noCPsEl
                RMtx(1,3*iCPs - 2) = dR(iCPs,1);
                RMtx(2,3*iCPs - 1) = dR(iCPs,1);
                RMtx(3,3*iCPs) = dR(iCPs,1);
            end
            
            %% 3iii.8. Compute the elementary area content at the Gauss Point
            elementaryAreaGP = detJCartesian2Param*detJParam2Integr*GW(iGP);
            
            %% 3iii.9. Compute and assemble the element load vector and the tangent matrix contribution at the Gauss Point
            F(EFT) = F(EFT) + RMtx'*tractionVct*elementaryAreaGP;
        end
    end
end

%% 4. Update by the existing load vector
if isvector(FOutdated)  
    F = F + FOutdated;   
end

%% 5. Appendix
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
