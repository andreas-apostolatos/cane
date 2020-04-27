function e = computeErrorStokesFlowInUnitSquareUnderBodyForces ...
    (p, q, Xi, Eta, CP, up, error)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the error corresponding to the Finite Element Method applied to
% the Stokes problem in a unit square domain under certain type of body
% forces.
%
%       Input :
%         p,q : Polynomial degrees
%      Xi,Eta : Knot vectors in xi,eta-direction
%          CP : Control point coordinates and weights
%          up : The nodal solution vector on the Control Points
%       error : Structure desicive on the error
%
%      Output :
%           e : The corresponding error
%
% Function layout :
%
% 0. Read input
%
% 1. Create the element nodal solution vectors
%
% 2. Get the quadrature for the integration of the error
%
% 3. Loop over all the elements (knot spans)
% ->
%    3i. Compute the determinant of the Jacobian to the transformation from the parameter space to the parent domain
%
%   3ii. Loop over all the Quadrature Points
%   ->
%        3ii.1. Compute the coordinates of the quadrature points at the NURBS domain via mapping
%
%        3ii.2. Issue the quadrature weights the compute the quadrature weight needed in the multiplication
%
%        3ii.3. Get the element nodal transport vector
%
%        3ii.4. Compute the selected component at the Gauss Point
%
%        3ii.5. Compute the Cartesian location of the Gauss Point
%
%        3ii.6. Compute the Jacobian to the transformation from the physical to the parameter space
%
%        3ii.7. Compute the analytical values at the Gauss Point location
%
%        3ii.8. Compute the error at the Gauss Point
%
%        3ii.9. Add the contribution from the Gauss Point
%   <-
% <-
% 
%% Function main body

%% 0. Read input

% Compute the control point coordinates for the deformed configuration
CPd = computeDisplacedControlPointsIGAIncompressibleFlow2D ...
    (CP, up, error.component);

% Number of knots in xi,eta-direction
numKnots_xi = length(Xi);
numKnots_eta = length(Eta);

% Number of Control Points in xi,eta-direction
numCPs_xi = length(CPd(:, 1, 1));
numCPs_eta = length(CPd(1, :, 1));

% Initialize output
e = 0;
eAnalytical = 0;

%% 1. Create the element nodal solution vectors

% Initialization over the array
uEl = zeros(numKnots_xi - p - 1, numKnots_eta - q - 1, 3*(p + 1)*(q + 1));

% Recursive assignment of the elements in the array
for j = (q+1):(numKnots_eta-q-1)
    for i = (p+1):(numKnots_xi-p-1)
        k=1; 
        for c = j-q-1:j-1 
            for b = i-p:i
                uEl(i,j,k)   = up(3*(c*numCPs_xi + b)-2);
                uEl(i,j,k + 1) = up(3*(c*numCPs_xi + b)-1);
                uEl(i,j,k + 2) = up(3*(c*numCPs_xi + b));
                k = k + 3;
            end
        end
    end
end

%% 2. Get the quadrature for the integration of the error

% Get the quadrature rule in u-direction
numGP_xi = error.int.uNGauss;
[GPxi, GW_xi] = gauss(numGP_xi);

% Get the quadrature rule in v-direction
numGP_eta = error.int.vNGauss;
[GP_eta, GW_eta] = gauss(numGP_eta);

%% 3. Loop over all the elements (knot spans)
for j = (q + 1):(numKnots_eta-q-1)   
    for i = (p+1):(numKnots_xi-p-1)
        % check if we are in a non-zero knot span
        if (Xi(i+1)~=Xi(i) && Eta(j+1)~=Eta(j))
            %% 3i. Compute the determinant of the Jacobian to the transformation from the parameter space to the parent domain
            detJxiu = (Xi(i+1)-Xi(i))*(Eta(j+1)-Eta(j))/4;
            
            %% 3ii. Loop over all the Quadrature Points
            for keta = 1:length(GP_eta)
                for kxi =1:length(GPxi)
                    %% 3ii.1. Compute the coordinates of the quadrature points at the NURBS domain via mapping
                    up = ( Xi(i+1)+Xi(i) + GPxi(kxi)*(Xi(i+1)-Xi(i)) )/2;
                    v = ( Eta(j+1)+Eta(j) + GP_eta(keta)*(Eta(j+1)-Eta(j)) )/2;
                    
                    %% 3ii.2. Issue the quadrature weights the compute the quadrature weight needed in the multiplication
                    gw = GW_xi(kxi)*GW_eta(keta);
                    
                    %% 3ii.3. Get the element nodal transport vector
                    uActual(:,1) = uEl(i,j,:);
                    
                    %% 3ii.4. Compute the selected component at the Gauss Point
                    upElement = computeNodalVectorIncompressibleFlow2D(i,p,up,Xi,j,q,v,Eta,CP,uActual);
                    if error.component==1||error.component==2||error.component==3
                        upFEM = upElement(error.component);
                    elseif error.component==4
                        upFEM = sqrt(upElement(1)^2 + upElement(2)^2);
                    end
                    
                    %% 3ii.5. Compute the Cartesian location of the Gauss Point
                    PCartesian = computePointCartesianCoordinatesOnBSplineSurface(p,i,up,Xi,q,j,v,Eta,CP);
                    x = PCartesian(1);
                    y = PCartesian(2);
                    
                    %% 3ii.6. Compute the Jacobian to the transformation from the physical to the parameter space
                    [~,dR] = computeNurbsBasisFunctionsAndFirstDerivatives2D(i,p,up,Xi,j,q,v,Eta,CP);
                    Jxxi = zeros(2,2);
                    k = 0;
                    for c = 0:q
                        for b = 0:p
                            % Update counter
                            k = k + 1;
                        
                            % Compute recursively the entries of the Jacobian
                            Jxxi(1,1) = Jxxi(1,1) + CP(i-p+b,j-q+c,1)*dR(k,1);
                            Jxxi(1,2) = Jxxi(1,2) + CP(i-p+b,j-q+c,2)*dR(k,1);
                            Jxxi(2,1) = Jxxi(2,1) + CP(i-p+b,j-q+c,1)*dR(k,2);
                            Jxxi(2,2) = Jxxi(2,2) + CP(i-p+b,j-q+c,2)*dR(k,2);
                        end
                    end
                    
                    %% 3ii.7. Compute the analytical values at the Gauss Point location
                    if error.component==1 || error.component==2 || error.component==4
                        % Compute the analytical value for the velocity field
                        if error.component==1
                            upAnalytical = x^2 * (1 - x)^2 * (2*y - 6*y^2 + 4*y^3);
                        elseif error.component==2
                            upAnalytical = -y^2 * (1 - y)^2 * (2*x - 6*x^2 + 4*x^3);
                        elseif error.component==1||error.component==4
                            ux = x^2 * (1 - x)^2 * (2*y - 6*y^2 + 4*y^3);
                            uy = -y^2 * (1 - y)^2 * (2*x - 6*x^2 + 4*x^3);
                            upAnalytical = sqrt(ux^2 + uy^2);
                        end
                    elseif error.component==3
                        upAnalytical = x*(1-x);
                    end
                    
                    %% 3ii.8. Compute the error at the Gauss Point
                    if strcmp(error.type, 'L2')
                        eGP = norm(upFEM - upAnalytical)^2*gw*detJxiu;
                        eAnalyticalGP = norm(upAnalytical)^2*gw*detJxiu;
                    end
                    
                    %% 3ii.9. Add the contribution from the Gauss Point
                    e = e + eGP;
                    eAnalytical = eAnalytical + eAnalyticalGP;
                end
            end
        end
    end
end

%% 4. Compute the relative error
e = e/eAnalytical;

end
