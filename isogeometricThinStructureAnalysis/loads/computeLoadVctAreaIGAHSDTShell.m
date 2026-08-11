function FGamma = computeLoadVctAreaIGAHSDTShell...
    (FGamma, BSplinePatch, xiLoadExtension, etaLoadExtension, loadAmplitude,...
    loadDirection, isConservative, t, int, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    [Your Name] (based on Andreas Apostolatos)
%
%% Function documentation
%
% Returns the consistent load vector for HSDT shell analysis with 5 DOF per
% control point corresponding to a distributed load over an area
%
%                Input :
%               FGamma : The existing load vector
%         BSplinePatch : The B-Spline patch structure
%        xiLoadExtension : Load extension in xi-direction 
%       etaLoadExtension : Load extension in eta-direction
%          loadAmplitude : The load amplitude
%          loadDirection : The direction of the applied load 
%                          ('x', 'y', 'z', or 'normal')
%        isConservative : Flag for conservative loading
%                     t : Time instance (for time-dependent loads)
%                   int : Integration structure
%                outMsg : Output message flag
%
%               Output :
%               FGamma : The updated load vector (5 DOF per control point)
%
% Function layout :
%
% 0. Read input
%
% 1. Choose integration scheme
%
% 2. Loop over elements in the load domain
%
% 3. Loop over Gauss points and compute load contributions
%
% 4. Assemble load vector
%
%% Function main body

%% 0. Read input

if strcmp(outMsg, 'outputEnabled')
    fprintf('_________________________________________________________________\n');
    fprintf('#################################################################\n');
    fprintf('Computation of the load vector for HSDT shell analysis has been \n');
    fprintf('initiated \n');
    fprintf('_________________________________________________________________\n\n');

    % start measuring computational time
    tic;
end

% Get patch properties
p = BSplinePatch.p;
q = BSplinePatch.q;
Xi = BSplinePatch.Xi;
Eta = BSplinePatch.Eta;  
CP = BSplinePatch.CP;
isNURBS = BSplinePatch.isNURBS;
DOFNumbering = BSplinePatch.DOFNumbering;

% Get load extension
xiSpan = xiLoadExtension;
etaSpan = etaLoadExtension;

% Number of knots and control points
numKnots_xi = length(Xi);
numKnots_eta = length(Eta);
numCPs_xi = length(CP(:, 1, 1));
numCPs_eta = length(CP(1, :, 1));

% Check if FGamma has correct size for HSDT (5 DOF per control point)
expectedSize = 5 * numCPs_xi * numCPs_eta;
if length(FGamma) ~= expectedSize
    error('Load vector size mismatch. Expected %d DOFs for HSDT (5 per CP), got %d', ...
          expectedSize, length(FGamma));
end

%% 1. Choose integration scheme

% Default integration
if strcmp(int.type, 'default')
    numGP_xi = p + 1;
    numGP_eta = q + 1;  
elseif strcmp(int.type, 'user')
    numGP_xi = int.xiNGPForLoad;
    numGP_eta = int.etaNGPForLoad;
end

% Get Gauss points and weights
[xiGP, xiGW] = getGaussPointsAndWeightsOverUnitDomain(numGP_xi);
[etaGP, etaGW] = getGaussPointsAndWeightsOverUnitDomain(numGP_eta);

%% 2. Loop over elements in the load domain

% Find element spans that intersect with load domain
xiSpanStart = max(1, find(Xi <= xiSpan(1), 1, 'last'));
xiSpanEnd = min(numKnots_xi - 1, find(Xi >= xiSpan(2), 1, 'first'));
etaSpanStart = max(1, find(Eta <= etaSpan(1), 1, 'last'));  
etaSpanEnd = min(numKnots_eta - 1, find(Eta >= etaSpan(2), 1, 'first'));

% Loop over element spans
for j = max(q + 1, etaSpanStart):min(numKnots_eta - q - 1, etaSpanEnd)
    for i = max(p + 1, xiSpanStart):min(numKnots_xi - p - 1, xiSpanEnd)
        if Xi(i + 1) ~= Xi(i) && Eta(j + 1) ~= Eta(j)
            
            % Element Jacobian determinant for parametric transformation
            detJxiu = (Xi(i + 1) - Xi(i)) * (Eta(j + 1) - Eta(j)) / 4;
            
            % Create element freedom table (5 DOF per control point)
            numDOFsEl = 5 * (p + 1) * (q + 1);
            EFT = zeros(1, numDOFsEl);
            k = 1;
            
            for cpj = j - q:j
                for cpi = i - p:i
                    EFT(k) = DOFNumbering(cpi, cpj, 1);     % u
                    EFT(k + 1) = DOFNumbering(cpi, cpj, 2); % v
                    EFT(k + 2) = DOFNumbering(cpi, cpj, 3); % w
                    EFT(k + 3) = DOFNumbering(cpi, cpj, 4); % θx
                    EFT(k + 4) = DOFNumbering(cpi, cpj, 5); % θy
                    k = k + 5;
                end
            end
            
            %% 3. Loop over Gauss points and compute load contributions
            for iGP_eta = 1:numGP_eta
                for iGP_xi = 1:numGP_xi
                    
                    % Gauss point coordinates in parametric space
                    xi = (Xi(i + 1) + Xi(i) + xiGP(iGP_xi) * (Xi(i + 1) - Xi(i))) / 2;
                    eta = (Eta(j + 1) + Eta(j) + etaGP(iGP_eta) * (Eta(j + 1) - Eta(j))) / 2;
                    
                    % Check if Gauss point is within load domain
                    if xi >= xiSpan(1) && xi <= xiSpan(2) && ...
                       eta >= etaSpan(1) && eta <= etaSpan(2)
                        
                        % Compute basis functions and derivatives
                        nDrvBasis = 1;
                        dR = computeIGABasisFunctionsAndDerivativesForSurface...
                            (i, p, xi, Xi, j, q, eta, Eta, CP, isNURBS, nDrvBasis);
                        
                        % Compute base vectors 
                        nDrvBaseVct = 0;
                        [A1, A2] = computeBaseVectorsAndDerivativesForBSplineSurface...
                            (i, p, j, q, CP, nDrvBaseVct, dR);
                        
                        % Surface normal and area element
                        A3Tilde = cross(A1(:, 1), A2(:, 1));
                        dA = norm(A3Tilde);
                        A3 = A3Tilde / dA;  % Normalized normal
                        
                        % Determine load direction vector
                        switch lower(loadDirection)
                            case 'x'
                                loadVec = [1; 0; 0];
                            case 'y'  
                                loadVec = [0; 1; 0];
                            case 'z'
                                loadVec = [0; 0; 1];
                            case 'normal'
                                loadVec = A3;  % Normal direction
                            otherwise
                                error('Unknown load direction: %s', loadDirection);
                        end
                        
                        % Element load vector (5 DOF per control point)
                        FeOnGP = zeros(numDOFsEl, 1);
                        
                        % Loop over element control points
                        for iCP = 1:(p + 1) * (q + 1)
                            % HSDT: only displacements receive direct loads
                            % Rotations are typically not directly loaded
                            
                            % Load contributions to displacements
                            FeOnGP(5*(iCP-1) + 1) = dR(iCP, 1) * loadAmplitude * loadVec(1); % u
                            FeOnGP(5*(iCP-1) + 2) = dR(iCP, 1) * loadAmplitude * loadVec(2); % v  
                            FeOnGP(5*(iCP-1) + 3) = dR(iCP, 1) * loadAmplitude * loadVec(3); % w
                            
                            % Rotational DOFs typically don't receive direct loads
                            % FeOnGP(5*(iCP-1) + 4) = 0;  % θx
                            % FeOnGP(5*(iCP-1) + 5) = 0;  % θy
                        end
                        
                        % Integration weight and Jacobian
                        weightGP = dA * detJxiu * xiGW(iGP_xi) * etaGW(iGP_eta);
                        
                        %% 4. Assemble load vector
                        FGamma(EFT) = FGamma(EFT) + FeOnGP * weightGP;
                        
                    end
                end
            end
        end
    end
end

%% 5. Appendix
if strcmp(outMsg, 'outputEnabled')
    computationalTime = toc;
    fprintf('Load vector computation for HSDT took %.2d seconds \n\n', computationalTime);
    fprintf('_____________________Load Vector Computation Ended_______________\n');
    fprintf('#################################################################\n\n\n');
end

end