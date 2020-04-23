function [F, tanMtx] = computeLoadVctAreaIGAThinStructure ...
    (FOutdated, BSplinePatch, xib, etab, FAmp, direction, isFollower, ...
    t, propInt, outMsg)
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
%        propInt : Structure containing information on the numerical
%                  quadrature
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
% 1. Get a quadrature rule
%
% 2. Find the start and the end knot span of the load application
%
% 3. Loop over the elements
% ->
%    3i. Compute the determinant of the Jacobian to the transformation from the parameter to the integration domain
%
%   3ii. Get the Element Freedom Table (EFT)
%
%  3iii. Loop over all the Gauss Points
%  ->
%        3iii.1. Compute the coordinates, the map and the Gauss Point location and weight for the fixed parametric coordinate
%
%        3iii.2. Compute the weight of the 2-dimensional quadrature
%
%        3iii.3. Compute the IGA basis functions and their first derivatives
%
%        3iii.4. Compute the base vectors
%
%        3iii.5. Compute the normal to the surface vector
%
%        3iii.6. Compute the traction vector
%
%        3iii.7. Compute the matrices containing the basis functions and their parametric derivatives
%
%        3iii.8. Formulate the cross product matrices g1x and g2x
%
%        3iii.9. Compute the elementary area content at the Gauss Point
%
%        3iii.10. Compute and assemble the element load vector and the tangent matrix contribution at the Gauss Point
%  <-
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
        fprintf('Update of the load vector corresponding to area\n');
    else
        fprintf('Computation of the load vector corresponding to area\n');
    end
    fprintf('load for the isogeometric Kirchhoff-Love shell problem\n');
    fprintf('has been initiated\n\n');
    if ~isa(FAmp, 'function_handle')
        fprintf('Constant area load is assumed with amplitude = %.2d\n', FAmp);
    else
        fprintf('Varying area load is assumed\n');
    end
    if strcmp(direction, 'x') || strcmp(direction, 'y') || strcmp(direction, 'z') ...
            || strcmp(direction, 'theta1') || strcmp(direction, 'theta2') ...
            || strcmp(direction, 'normal')
        fprintf('Load direction is chosen as %s', direction);
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
numCPs_xi = length(CP(:, 1, 1));
numCPs_eta = length(CP(1, :, 1));
if isFollower
    CPd = BSplinePatch.CPd;
end
isNURBS = BSplinePatch.isNURBS;

% Number of Control Points affecting each element
noCPsEl = (p + 1)*(q + 1);

% Number of DOFs per Control Point
numDOFsCP = 3;

% Number of DOFs per element
numDOFsEl = numDOFsCP*noCPsEl;

% Numner of derivatives for the B-Spline basis functions
numDrvsBasisFct = 1;

% Number of derivatives for the base vectors
numDrvsBaseVcts = 0;

% Number of DOFs for the patch
if isfield(BSplinePatch, 'noDOFs')
    noDOFs = BSplinePatch.noDOFs;
else
    noDOFs = numDOFsCP*numCPs_xi*numCPs_eta;
end

% Initialize auxiliary arrays
RMtx = zeros(numDOFsCP, numDOFsEl);
if isFollower
    dRdXiMtx = zeros(numDOFsCP, numDOFsEl);
    dRdEtaMtx = zeros(numDOFsCP, numDOFsEl);
end

% Initialize output arrays
F = zeros(noDOFs, 1);
if isFollower
    tanMtx = zeros(noDOFs);
else
    tanMtx = 'undefined';
end

%% 1. Get a quadrature rule
if isstruct(propInt)
    if ~isfield(propInt, 'type')
        error('propInt must define a type')
    end
    if strcmp(propInt.type,'default')
        numGP_xi = ceil(2*p - 1);
        numGP_eta = ceil(2*q - 1);
    elseif strcmp(propInt.type,'user')
        numGP_xi = propInt.xiNGPForLoad;
        numGP_eta = propInt.etaNGPForLoad;
    else
        error('Select a valid integration type int.type')
    end
else
    error('propInt must be a structure');
end
[GP_xi,GW_xi] = getGaussPointsAndWeightsOverUnitDomain(numGP_xi);
[GP_eta,GW_eta] = getGaussPointsAndWeightsOverUnitDomain(numGP_eta);

%% 2. Find the start and the end knot span of the load application

% On the xi-direction
xiSpanStart = findKnotSpan(xib(1), Xi, numCPs_xi);
xiSpanEnd = findKnotSpan(xib(2), Xi, numCPs_xi);

% On the eta-direction
etaSpanStart = findKnotSpan(etab(1), Eta, numCPs_eta);
etaSpanEnd = findKnotSpan(etab(2), Eta, numCPs_eta);

%% 3. Loop over the elements
for iEtaSpan = etaSpanStart:etaSpanEnd
    for iXiSpan = xiSpanStart:xiSpanEnd
        % check if element is greater than zero
        if Xi(iXiSpan + 1) ~= Xi(iXiSpan) && Eta(iEtaSpan + 1) ~= Eta(iEtaSpan)
            %% 3i. Compute the determinant of the Jacobian to the transformation from the parameter to the integration domain
            detJParam2Integr = (Xi(iXiSpan + 1) - Xi(iXiSpan))*(Eta(iEtaSpan + 1) - Eta(iEtaSpan))/4;
            
            %% 3ii. Get the Element Freedom Table (EFT)
            idElmt = BSplinePatch.knotSpan2ElmntNo(iXiSpan, iEtaSpan);
            EFT = BSplinePatch.EFT(:, idElmt);
            
            %% 3iii. Loop over all the Gauss Points
            for iGPEta = 1:numGP_eta
                for iGPXi= 1:numGP_xi
                    %% 3iii.1. Compute the coordinates, the map and the Gauss Point location and weight for the fixed parametric coordinate
                    xi = (Xi(iXiSpan + 1) + Xi(iXiSpan) + GP_xi(iGPXi)*(Xi(iXiSpan + 1) - Xi(iXiSpan)))/2;
                    eta = (Eta(iEtaSpan + 1) + Eta(iEtaSpan) + GP_eta(iGPEta)*(Eta(iEtaSpan + 1) - Eta(iEtaSpan)))/2;
                    
                    %% 3iii.2. Compute the weight of the 2-dimensional quadrature
                    GW = GW_xi(iGPXi)*GW_eta(iGPEta);
                    
                    %% 3iii.3. Compute the IGA basis functions and their first derivatives
                    dR = computeIGABasisFunctionsAndDerivativesForSurface ...
                        (iXiSpan, p, xi, Xi, iEtaSpan, q, eta, Eta, CP, ...
                        isNURBS, numDrvsBasisFct);
                    
                    %% 3iii.4. Compute the base vectors
                    if isFollower
                        [a1, a2] = computeBaseVectorsAndDerivativesForBSplineSurface ...
                            (iXiSpan, p, iEtaSpan, q, CPd, numDrvsBaseVcts, dR);
                    else
                        [A1, A2] = computeBaseVectorsAndDerivativesForBSplineSurface ...
                            (iXiSpan, p, iEtaSpan, q, CP, numDrvsBaseVcts, dR);
                    end
                    
                    %% 3iii.5. Compute the normal to the surface vector
                    if isFollower
                        nTilde = cross(a1, a2);
                        detJCartesian2Param = norm(nTilde);
                        n = nTilde/detJCartesian2Param;
                    else
                        NTilde = cross(A1, A2);
                        detJCartesian2Param = norm(NTilde);
                        N = NTilde/detJCartesian2Param;
                    end
                    
                    
                    %% 3iii.6. Compute the traction vector
                    
                    % Compute the amplitude of the force if it is spatially 
                    % dependent
                    if ~isnumeric(FAmp)
                        X = computeCartesianCoordinatesOfAPointOnBSplineSurface ...
                            (iXiSpan, p, xi, Xi, iEtaSpan, q, eta, Eta, CP, dR(:, 1));
                        if  ischar(FAmp) && ~isa(FAmp,'function_handle')
                            FAmp = str2func(FAmp);
                        end
                        FAmplitude = FAmp(X(1, 1), X(2, 1),X(3, 1), t);
                    else
                        FAmplitude = FAmp;
                    end
                    
                    % Decide upon the direction
                    if strcmp(direction, 'x')
                        dirVct = [1
                                  0 
                                  0];
                    elseif strcmp(direction, 'y')
                        dirVct = [0 
                                  1 
                                  0];
                    elseif strcmp(direction, 'z')
                        dirVct = [0 
                                  0 
                                  1];
                    elseif strcmp(direction, 'theta1')
                        if isFollower
                            dirVct = a1/norm(a1);
                        else
                            dirVct = A1/norm(A1);
                        end
                    elseif strcmp(direction, 'theta2')
                        if isFollower
                            dirVct = a2/norm(a2);
                        else
                            dirVct = A2/norm(A2);
                        end
                    elseif strcmp(direction, 'normal')
                        if isFollower
                            dirVct = n;
                        else
                            dirVct = N;
                        end
                    end
                    
                    % Compute the tractionvector
                    tractionVct = FAmplitude*dirVct;
                    
                    %% 3iii.7. Compute the matrices containing the basis functions and their parametric derivatives
                    for iCPs = 1:noCPsEl
                        RMtx(1, 3*iCPs - 2) = dR(iCPs, 1);
                        RMtx(2, 3*iCPs - 1) = dR(iCPs, 1);
                        RMtx(3, 3*iCPs) = dR(iCPs, 1);
                        if isFollower && strcmp(direction, 'normal')
                            dRdXiMtx(1, 3*iCPs - 2) = dR(iCPs, 2);
                            dRdXiMtx(2, 3*iCPs - 1) = dR(iCPs, 2);
                            dRdXiMtx(3, 3*iCPs) = dR(iCPs, 2);
                            dRdEtaMtx(1, 3*iCPs - 2) = dR(iCPs, 3);
                            dRdEtaMtx(2, 3*iCPs - 1) = dR(iCPs, 3);
                            dRdEtaMtx(3, 3*iCPs) = dR(iCPs, 3);
                        end
                    end
                    
                    %% 3iii.8. Formulate the cross product matrices g1x and g2x
                    if isFollower && strcmp(direction, 'normal')
                        a1x = [0         -a1(3, 1) a1(2, 1)
                               a1(3, 1)  0        -a1(1, 1)
                               -a1(2, 1) a1(1, 1) 0];
                        a2x = [0         -a2(3, 1) a2(2, 1)
                               a2(3, 1)  0        -a2(1, 1)
                               -a2(2, 1) a2(1, 1)  0];
                    end
                    
                    %% 3iii.9. Compute the elementary area content at the Gauss Point
                    elementaryAreaGP = detJCartesian2Param*detJParam2Integr*GW;
                    
                    %% 3iii.10. Compute and assemble the element load vector and the tangent matrix contribution at the Gauss Point
                    F(EFT) = F(EFT) + (RMtx'*tractionVct)*elementaryAreaGP;
%                     F(EFT) = F(EFT) + ((RMtx'*tractionVct) + FAmp*(-a2x*dRdXiMtx + a1x*dRdEtaMtx)'*RMtx*dHat(EFT))*elementaryAreaGP;

                   if isFollower && strcmp(direction,'normal')
                       tanMtx(EFT,EFT) = tanMtx(EFT,EFT) - ...
                           FAmplitude*(RMtx'*(-a2x*dRdXiMtx + a1x*dRdEtaMtx))*detJParam2Integr*GW;
                   end
                end
            end
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
