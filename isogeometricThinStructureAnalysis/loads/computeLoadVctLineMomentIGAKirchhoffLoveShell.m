function [F, tanMtx] = computeLoadVctLineMomentIGAKirchhoffLoveShell ...
    (FOutdated, BSplinePatch, xib, etab, MAmp, direction, isFollower, ...
    t, int, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the consistent nodal forces corresponding to a bending monent
% distributed over the edge of the shell. The direction of the externally
% applied moment can be in the normal to the shell boundary direction (i.e. 
% bending moment) or in the tangent to the shell boundary direction (i.e 
% twisting momemt)
%
%          Input :
%      FOutdated : Existing load vector
%   BSplinePatch : Array containing information on the B-Spline patch
%                     .p,.q : Polynomial orders in xi- and -eta direction
%                  .Xi,.Eta : Polynomial orders in xi-eta direction
%                       .CP : Set of Control Point coordinates and weights
%                  .isNURBS : Flag on whether the basis is a NURBS or a
%                             B-Spline
%                   .noDOFs : Total Number of DOFs for the B-Spline patch
%                             including possible Lagrange Multipliers
%       xib,etab : Load extension (e.g. xib = [0 1], etab = 1)
%           MAmp : constant line moment or handle to moment load function
%      direction : direction of the moment MAmp:
%                       'bending' : The applied moment is a bending moment
%                      'twisting' : The applied moment is a twisting moment
%     isFollower : Flag on whether the applied load is a follower load
%              t : The time instance where the load vector is computed
%            int : structure deciding upon the integration scheme
%         outMsg : Whether or not to output message on refinement progress
%                   'outputEnabled' : enables output information
%
%         Output :
%              F : The updated force vector
%         tanMtx : The tangent matrix from the application of a follower
%                  load
%
% Function layout:
%
% 0. Read input
%
% 1. Get the running and the fixed parameter on the patch interface and the coupling region
%
% 2. Issue Gauss Point coordinates and weights
%
% 3. Loop over all the elements on the boundary
% ->
%    3i. Compute the map from the parameter to the integration space
%
%   3ii. Get the Element Freedom Table (EFT)
%
%  3iii. Loop over all Gauss points
%  ->
%        3iii.1. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
%
%        3iii.2. Compute the NURBS basis functions
%
%        3iii.3. Compute the covariant base vectors
%
%        3iii.4. Compute the surface normal vector
%
%        3iii.5. Compute the derivatives of the surface normal vector
%
%        3iii.6. Compute the normal to the boundary vector
%
%        3iii.7. Compute the covariant metric coefficients
%
%        3iii.8. Compute the contravariant base vectors
%
%        3iii.9. Compute the basis functions matrix and their derivatives
%
%        3iii.10. Transform the normal and the tangent vector to the contravariant basis
%
%        3iii.11. Compute the curvature coefficients
%
%        3iii.12. Compute the B-operator matrix for the rotations
%
%        3iii.13. Compute the Jacobians for the transformations to the parent domain
%
%        3iii.14. On the applied moment
%
%        3iii.15. Compute and assemble the local force vector at the Gauss point to the global load vector
%  <-
% <-
%
% 4. Update the load vector
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
    fprintf('boundary moment for the isogeometric Kirchhoff-Love shell\n');
    fprintf('problem has been initiated\n\n');
    if isnumeric(MAmp)
        fprintf('Constant boundary load is assumed with amplitude = %.2d\n',MAmp);
    else
        fprintf('Varying boundary load is assumed\n');
    end
    if strcmp(direction,'bending')
        fprintf('Application of a boundary bending moment\n');
    elseif strcmp(direction,'twisting')
        fprintf('Application of a boundary twisting moment\n');
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
if ~isFollower
    CP = BSplinePatch.CP;
else
    CP = BSplinePatch.CPd;
end
isNURBS = BSplinePatch.isNURBS;
if isfield(BSplinePatch,'noDOFs')
    noDOFs = BSplinePatch.noDOFs;
else
    noDOFs = 3*length(CP(:,1,1))*length(CP(1,:,1));
end

% Number of control points in u,v-direction
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));

% Number of local DOFs
noCPsEl = (p + 1)*(q + 1);
noDOFsLoc = 3*noCPsEl;

% Initialize auxiliary arrays
RMtx = zeros(3,noDOFsLoc);
dRdXiMtx = zeros(3,noDOFsLoc);
dRdEtaMtx = zeros(3,noDOFsLoc);

% Initialize output array
F = zeros(noDOFs,1);
if isFollower
    tanMtx = zeros(noDOFs);
else
    tanMtx = 'undefined';
end

%% 1. Get the running and the fixed parameter on the patch interface and the coupling region
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

%% 2. Issue Gauss Point coordinates and weights
if strcmp(int.type,'default')
    if isscalar(etab)
        noGPs = p + 1;
    elseif isscalar(xib)
        noGPs = q + 1;
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

%% 3. Loop over all the elements on the boundary
for iKnotSpan = thetaStartKnotSpan:thetaEndKnotSpan
    % check if we are in a non-zero knot span
    if thetaKnotVct(iKnotSpan + 1) ~= thetaKnotVct(iKnotSpan)
        %% 3i. Compute the map from the parameter to the integration space
        detJxizeta = (thetaKnotVct(iKnotSpan + 1) - thetaKnotVct(iKnotSpan))/2;
        
        %% 3ii. Get the Element Freedom Table (EFT)
        if isOnXi
            xiKnotSpan = iKnotSpan;
        else
            etaKnotSpan = iKnotSpan;
        end
        noElmnt = BSplinePatch.knotSpan2ElmntNo(xiKnotSpan,etaKnotSpan);
        EFT = BSplinePatch.EFT(:,noElmnt);
        
        %% 3iii. Loop over all Gauss points
        for iGP = 1:noGPs
            %% 3iii.1. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
            theta = ((1 - GP(iGP))*thetaKnotVct(iKnotSpan) + (1 + GP(iGP))*thetaKnotVct(iKnotSpan + 1))/2;
            
            %% 3iii.2. Compute the NURBS basis functions
            if isOnXi
                xi = theta;
                xiKnotSpan = findKnotSpan(xi,Xi,nxi);
            else
                eta = theta;
                etaKnotSpan = findKnotSpan(eta,Eta,neta);
            end
            noDrv = 2;
            dR = computeIGABasisFunctionsAndDerivativesForSurface...
                (xiKnotSpan,p,xi,Xi,etaKnotSpan,q,eta,Eta,CP,isNURBS,noDrv);
            
            %% 3iii.3. Compute the covariant base vectors
            noDrv = 1;
            [dA1,dA2] = computeBaseVectorsAndDerivativesForBSplineSurface...
                (xiKnotSpan,p,etaKnotSpan,q,CP,noDrv,dR);
            
            %% 3iii.4. Compute the surface normal vector
            A3Tilde = cross(dA1(:,1),dA2(:,1));
            A3 = A3Tilde/norm(A3Tilde);
            
            %% 3iii.5. Compute the derivatives of the surface normal vector
            dA3 = computeParametricDrvsSurfaceNormalOnBSplineSurface...
                ([dA1(:,1) dA2(:,1)],[dA1(:,2) dA2(:,2) dA1(:,3)],...
                A3,norm(A3Tilde));
            
            %% 3iii.6. Compute the normal to the boundary vector
            [nCartesian,tCartesian] = computeNormalAndTangentVectorsToBSplineBoundary...
                (xi,Xi,eta,Eta,dA1(:,1),dA2(:,1),A3,isOnXi);
            
            %% 3iii.7. Compute the covariant metric coefficients
            AabCov = [dA1(:,1) dA2(:,1)]'*[dA1(:,1) dA2(:,1)];

            %% 3iii.8. Compute the contravariant base vectors
            AContravariant = (AabCov\[dA1(:,1) dA2(:,1)]')';
            
            %% 3iii.9. Compute the basis functions matrix and their derivatives
            for iCPs = 1:noCPsEl
                RMtx(1,3*iCPs-2) = dR(iCPs,1);
                RMtx(2,3*iCPs-1) = dR(iCPs,1);
                RMtx(3,3*iCPs) = dR(iCPs,1);
                dRdXiMtx(1,3*iCPs-2) = dR(iCPs,2);
                dRdXiMtx(2,3*iCPs-1) = dR(iCPs,2);
                dRdXiMtx(3,3*iCPs) = dR(iCPs,2);
                dRdEtaMtx(1,3*iCPs-2) = dR(iCPs,4);
                dRdEtaMtx(2,3*iCPs-1) = dR(iCPs,4);
                dRdEtaMtx(3,3*iCPs) = dR(iCPs,4);
            end

            %% 3iii.10. Transform the normal and the tangent vector to the contravariant basis
            nContravariant = AContravariant'*nCartesian;
            tContravariant = AContravariant'*tCartesian;
            
            %% 3iii.11. Compute the curvature coefficients
            BV = [dA1(:,2) dA2(:,2) dA1(:,3)]'*(A3Tilde/norm(A3Tilde));

            %% 3iii.12. Compute the B-operator matrix for the rotations
            [Bt,Bn,~,~] = computeBOperatorMatrix4RotationsIGAKirchhoffLoveShell(RMtx,dRdXiMtx,dRdEtaMtx,...
                A3,dA3,AContravariant,BV,nContravariant,tContravariant);

            %% 3iii.13. Compute the Jacobians for the transformations to the parent domain
            if isOnXi
                detJxxi = norm(dA1(:,1));
            else
                detJxxi = norm(dA2(:,1));
            end
            
            %% 3iii.14. On the applied moment
            if isnumeric(MAmp)
                MAmplitude = MAmp;
            else
                MAmplitude = MAmp(x,y,z,t);
            end

            %% 3iii.15. Compute and assemble the local force vector at the Gauss point to the global load vector
            if strcmp(direction,'bending')
                F(EFT) = F(EFT) + Bt'*MAmplitude*detJxxi*detJxizeta*GW(iGP);
            elseif strcmp(direction,'twisting')
                F(EFT) = F(EFT) + Bn'*MAmplitude*detJxxi*detJxizeta*GW(iGP);
            end
        end
    end
end

%% 4. Update the load vector
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
