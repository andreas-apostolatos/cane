function [Fl, tanMtxLoadPatch] = ...
    computeLoadVctLineIGAIncompressibleNavierStokesFlow ...
    (FlOld, BSplinePatch, xib, etab, hFlux, direction, isFollower, ...
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
% Returns the consistent force vector for the application of a flux over a
% portion of the computational domain's boundary (Neumann boundary) defined
% by parameters ub, vb and which is of magnitude hFlux for a 2D
% incompressile flow problem.
% 
%        Input :
%        FlOld : The existing force vector
% BSplinePatch : The polynomial degrees, the knot vectors and the Control
%                Point coordinates/weights corresponding to the underlying
%                B-Spline patch
%     xib,etab : Flux application extension (e.g. ub=[0 1], vb=1)
%        hFlux : constant line flux or handle to load function [N/m]
%    direction : direction of h:
%                1=x, 2=y for h(s)   -> integration over the boundary space
%                3=x, 4=y for h(x,y) -> integration over the Cartesian 
%                                       space
%                5=parallel, 6=+90  to the edge
%                xiEta: dir to integrate: 1=xi, 2=eta determined from xib, 
%                       etab
%                  xys: 1=integrate ds   2=integrate dx,dy.  determined 
%                       from dir
%   isFollower : Flag on whether the load is a follower or not
%            t : The current time
%          int : Structure responsible for the integration
%               default : Automatic choice of the Gauss Points fulfiling
%                         stability
%                manual : Manual choice of the Gauss Points for improving
%                         the accuracy
%       outMsg : Whether or not to output message
%                'outputEnabled' : enables output information
%
%       Output :
%           Fl : The updated force vector
%
% Function layout :
%
% 0. Read input
%
% 1. On the numerical integration
%
% 2. Find the parametrization of the line integral
%
% 3. Loop over all elements of the load application
% ->
%    3i. Loop over all Gauss points
%    ->
%        3i.1. Compute the coordinate, the map and the Gauss Point location and weight for the fixed parametric coordinate
%
%        3i.2. Compute the IGA basis functions and their first derivatives
%
%        3i.3. Compute the the Jacobian of the transformation from the surface physical to the parameter space and the physical location of the load application
%
%        3i.4. Compute the load amplitude on the Gauss Point is the load is not constant
%
%        3i.5. Compute the product R(xi,eta)*ds in the parameter space or R(xi,eta)*dX in the physical space at the Gauss Point
%
%        3i.6. Compute the determinant of the Jacobian from the parameter to the integration space and the Gauss Point weight
%
%        3i.7. Compute the element load vector
%
%        3i.8. Assemble the element load vector to the global load vector
%    <-
% <-
%
% 4. Re-assemble the load vector with respect to the global numbering
%
% 5. Update the existing load vector
%
% 6. Appendix
%
%% Function main body
if strcmp(outMsg, 'outputEnabled')
    fprintf('_______________________________________________________\n');
    fprintf('#######################################################\n');
    if isvector(FlOld)
        fprintf('Update of the load vector corresponding to line\n');
    else
        fprintf('Computation of the load vector corresponding to line\n');
    end
    fprintf('boundary fluxes for the incompressible Navier-Stokes\n');
    fprintf('flow problem has been initiated\n\n');
    if isnumeric(MAmp)
        fprintf('Constant boundary flux is assumed with amplitude = %.2d\n', MAmp);
    else
        fprintf('Varying boundary flux is assumed\n');
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

% Get B-Spline patch parameters
p = BSplinePatch.p;
q = BSplinePatch.q;
Xi = BSplinePatch.Xi;
Eta = BSplinePatch.Eta;
CP = BSplinePatch.CP;
isNURBS = BSplinePatch.isNURBS;

% Length of the knot vectors
numKnots_xi = length(Xi);
numKnots_eta = length(Eta);

% Number of control points in u,v-direction
numCPs_xi = length(CP(:, 1, 1));
numCPs_neta = length(CP(1, :, 1));

% if the fload parameter is numeric then assign directly its value
if isnumeric(hFlux)==1
    h = hFlux;
end

% Number of DOFs
numDOFs = 2*numCPs_xi*numCPs_neta;

% Initialize auxiliary array
Rds = zeros(p + 1, q + 1);

% Initialize global force vector
F = zeros(numCPs_xi, numCPs_neta, 3);   

% Initialize element load vector
Fel = zeros(p + 1, q + 1, 3);

% Initialize output arrays
Fl = zeros(numDOFs, 1);
tanMtxLoadPatch = 'undefined';

%% 1. On the numerical integration

% Initialize structure
nGP = zeros(2, 1);

% Issue the Gauss points for the selected integration scheme:
if strcmp(int.type, 'default')
   % Default scheme is the full gaussian quadrature element-wise (FGI)
   nGP(1) = ceil((p + 1)/2) ;
   nGP(2) = ceil((q + 1)/2) ;
elseif strcmp(int.type, 'manual')
    % Manual choice of the gauss points
    nGP(1) = int.xiNGPForLoad;
    nGP(2) = int.etaNGPForLoad;
end

% Issue the Gauss points for the numerical integration
[xiGP,xiGW] = getGaussPointsAndWeightsOverUnitDomain(nGP(1));
[etaGP,etaGW] = getGaussPointsAndWeightsOverUnitDomain(nGP(2));

%% 2. Find the parametrization of the line integral

% For the xi-direction :
% ______________________

% Find the span in xi-direction where to apply the load
i1 = findKnotSpan(xib(1), Xi, numCPs_xi);

% Check if the load extension in the xi-direction has zero measure
if xib(1) == xib(2)
    % If its a scalar one Gauss point sufficies
    xiNGP = 1;   
    xigw = 1;
    
    % Get coordinate for the integration
    xi = xib(1);    
    
    % The start and the end span are identical
    i2 = i1;
    
    % The map from the parameter to the integration space is one
    xiMap = 1;     
    
    % Integration takes place along eta parametric coordinate
    xiEta = 2;
else
    % If xib is a vector then find the knot span of the end coordinate in
    % xib
    i2 = findKnotSpan(xib(2), Xi, numCPs_xi);
    
    % If the end span of the integration domain is not the end of the knot
    % span decrease the knot span index by 1 (I have no clue why he 
    % implemented that in this way)
    if xib(2) ~= Xi(numKnots_xi)
        i2 = i2 - 1;    
    end
    
    % Get the number of Gauss points for the integration along
    % xi-coordinate
    xiNGP = nGP(1, 1);
end

% For the eta-direction :
% _______________________

% Find the span in v-direction where to apply the load
j1 = findKnotSpan(etab(1), Eta, numCPs_neta);

% Check if the load extension in the eta-direction has zero measure
if etab(1) == etab(2)
    % If its a scalar one Gauss point sufficies
    etaNGP = 1;
    etagw = 1;
    
    % Get coordinate for the integration
    eta = etab(1);
    
    % The start and the end span are identical
    j2 = j1;
    
    % The map from the parameter to the integration space is one
    etaMap = 1;     
    
    % Integration takes place along xi parametric coordinate
    xiEta = 1;
else
    % If etab is a vector then find the knot span of the end coordinate in
    % etab
    j2 = findKnotSpan(etab(2), Eta, numCPs_neta);
    
    % If the end span of the integration domain is not the end of the knot
    % span decrease the knot span index by 1 (I have no clue why he 
    % implemented that in this way)
    if etab(2) ~= Eta(numKnots_eta)
        j2 = j2 - 1;    
    end
    
    % Get the number of Gauss points for the integration along
    % eta-coordinate
    etaNGP = nGP(2, 1);
end

% Decide on which coordinate system to perform the integration
if  direction == 1 || direction == 2 || direction == 3|| direction == 7 || ...
        direction == 8 || direction == 9
    % Integrate over the curvilinear parametric line
    xys = 1;
elseif direction == 4 || direction == 5 || direction == 6
    % Integrate over the cartesian coordinate system xy
    xys = 2;
end

%% 3. Loop over all elements of the load application
for j = j1:j2
    for i = i1:i2
        % check if we are in a non-zero knot span
        if Xi(i+1)~=Xi(i) && Eta(j+1)~=Eta(j)
            %% 3i. Loop over all Gauss points
            for cEta = 1:etaNGP
                for cXi = 1:xiNGP
                    %% 3i.1. Compute the coordinate, the map and the Gauss Point location and weight for the fixed parametric coordinate

                    % For the xi-direction :
                    % ______________________

                    if ~isscalar(xib)
                        % map the quadrature point into the knot span in
                        % xi-direction
                        xi = (Xi(i+1) + Xi(i) + xiGP(cXi)*(Xi(i + 1) - Xi(i)))/2;

                        % Compute the respective Jacobian determinant
                        xiMap = (Xi(i + 1) - Xi(i))/2;

                        % Issue quadrature weight in u-direction
                        xigw = xiGW(cXi);
                    end

                    % For the eta-direction :
                    % _______________________

                    if ~isscalar(etab)
                        % map the quadrature point into the knot span in
                        % eta-direction
                        eta = (Eta(j + 1) + Eta(j) + etaGP(cEta)*(Eta(j + 1) - Eta(j)))/2;

                        % compute the respective Jacobian determinant
                        etaMap = (Eta(j + 1) - Eta(j))/2;

                        % issue quadrature weight in v-direction
                        etagw = etaGW(cEta);
                    end

                    %% 3i.2. Compute the IGA basis functions and their first derivatives
                    numDrv = 1;
                    dR = computeIGABasisFunctionsAndDerivativesForSurface ...
                        (i, p, xi, Xi, j, q, eta, Eta, CP, isNURBS, numDrv);

                    %% 3i.3. Compute the the Jacobian of the transformation from the surface physical to the parameter space and the physical location of the load application
                    
                    % Initialize the Jacobian matrix
                    J = zeros(2, 2);
                    
                    % Initialize the Cartesian coordinates of the load
                    % application point if the load is not constant but
                    % varying
                    if isnumeric(hFlux) == 0
                        x = 0;
                        y = 0;
                    end
                    
                    % Compute the entries of the Jacobian matrix
                    % iteratively
                    counter = 0;
                    for b = 0:q
                        for a = 0:p
                            % Update counter
                            counter = counter + 1;
                            
                            % Compute the entries of the Jacobian matrix
                            J(1, 1) = J(1, 1) + CP(i - p + a, j - q + b, 1)*dR(counter, 2);
                            J(1, 2) = J(1, 2) + CP(i - p + a, j - q + b, 2)*dR(counter, 2);
                            J(2, 1) = J(2, 1) + CP(i - p + a, j - q + b, 1)*dR(counter, 3);
                            J(2, 2) = J(2, 2) + CP(i - p + a, j - q + b, 2)*dR(counter, 3);
                            
                            % Compute the physical location of the load
                            % application point if the load is not constant
                            if isnumeric(hFlux) == 0
                                x = dR(counter, 1)*CP(i - p + a, j - q + b, 1) + x;
                                y = dR(counter, 1)*CP(i - p + a, j - q + b, 2) + y;
                            end
                        end
                    end
                    
                    %% 3i.4. Compute the load amplitude on the Gauss Point is the load is not constant
                    if isnumeric(hFlux) == 0
                        h = hFlux(p, i, xi, Xi, q, j, eta, Eta, CP, t);  
                    end
                    
                    %% 3i.5. Compute the product R(xi,eta)*ds in the parameter space or R(xi,eta)*dX in the physical space at the Gauss Point
                    
                    % Initialize counter
                    counter = 0;
                    
                    % If the integration is performed over the parameter
                    % space compute R(xi,eta)*ds
                    if xys == 1
                        % If we integrate along the parametric line the 
                        % determinant of the Jacobian is the velocity of 
                        % the parametrization in the NURBS space
                        dsdxiEta = sqrt(J(xiEta, 1)^2 + J(xiEta, 2)^2);

                        for b = 0:q
                            for a = 0:p
                                % Update counter
                                counter = counter + 1;
                                
                                % Compute the product of the basis function
                                % with the determinant of the Jacobian
                                Rds(a + 1, b + 1) = dR(counter, 1)*dsdxiEta;
                            end
                        end
                    % Else if the integration is performed over the parameter
                    % space compute R(xi,eta)*dX
                    elseif xys == 2
                        for b = 0:q
                            for a = 0:p
                                % Update counter
                                counter = counter + 1;
                                Rds(a + 1, b + 1, 1) = dR(counter, 1)*J(xiEta, 1);
                                Rds(a + 1, b + 1, 2) = dR(counter, 1)*J(xiEta, 2);
                            end
                        end
                    end
                    
                    %% 3i.6. Compute the determinant of the Jacobian from the parameter to the integration space and the Gauss Point weight
                    
                    % Compute the determinant of the Jacobian from the 
                    % physical to the parameter space
                    detJParam2Integr = xiMap*etaMap;
                    
                    % Compute the product of the quadrature weights 
                    GW = xigw*etagw;
                    
                    %% 3i.7. Compute the element load vector
                    if direction == 1 || direction == 2 || direction == 3
                        Fel(:, :, direction) = h*Rds(:, :)*GW*detJParam2Integr;
                    elseif direction == 4
                        Fel(:, :, 1) = h*abs(Rds(:, :, 2))*GW*detJParam2Integr;
                    elseif direction == 5
                        Fel(:, :, 2) = h*abs(Rds(:, :, 1))*GW*detJParam2Integr;
                    elseif direction == 6
                        Fel(:, :, 1) = h*Rds(:, :, 1)*GW*detJParam2Integr;
                        Fel(:, :, 2) = h*Rds(:, :, 2)*GW*detJParam2Integr;
                    end
                    
                    %% 3i.8. Assemble the element load vector to the global load vector
                    F(i - p:i, j - q:j, :) = Fel(:, :, :) + F(i - p:i, j - q:j, :);
                end  
            end
        end
    end
end

%% 4. Re-assemble the load vector with respect to the global numbering
counter = 1;
for j = 1:length(F(1, :, 1))
    for i = 1:length(F(:, 1, 1))
        % Assemble the x-coordinates of the load vector
        Fl(counter, 1) = F(i, j, 1);
        
        % Assemble the y-coordinates of the load vector
        Fl(counter + 1, 1) = F(i, j, 2);
        
        % Assemble the z-coordinates of the load vector
        Fl(counter + 2, 1) = F(i, j, 3);
        
        % Update counter
        counter = counter + 3;
    end
end

%% 5. Update the existing load vector
if isvector(FlOld)  
    Fl = Fl + FlOld;   
end

%% 6. Appendix
if strcmp(outMsg, 'outputEnabled')
    computationalTime = toc;
    if isvector(FlOld)
        fprintf('Load vector update took %.2d seconds \n\n', computationalTime);
        fprintf('________________Load Vector Update Ended_______________\n');
    else
        fprintf('Load vector computation took %.2d seconds \n\n', computationalTime);
        fprintf('____________Load Vector Computation Ended______________\n');
    end
    fprintf('#######################################################\n\n\n');
end

end
