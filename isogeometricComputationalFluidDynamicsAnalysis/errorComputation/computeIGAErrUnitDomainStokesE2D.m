function [errL2Velocity, errL2Pressure, minElASize] = ...
    computeIGAErrUnitDomainStokesE2D ...
    (p, Xi, q, Eta, CP, isNURBS, up, propInt)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the L2-error of the velocity field and the L2-error of the
% pressure field for the benchmark problem of the Stokes flow in unit
% square.
%
%         Input :
%           p,q : The polynomial degrees in xi-,eta-direction
%        Xi,Eta : The knot vectors in xi-,eta-direction
%            CP : The set of Control Point coordinates and weights
%       isNURBS : Flag on whether the basis is a NURBS or a B-Spline
%       propInt : On the spatial integration
%                       .type : 'default' or 'manual'
%                      .xiNGP : No. of GPs along xi-direction for stiffness 
%                               entries
%                     .etaNGP : No. of GPs along eta-direction for 
%                               stiffness entries
%               .xiNGPForLoad : No. of GPs along xi-direction for load entries
%              .etaNGPForLoad : No. of GPs along eta-direction for load entries
%                 .nGPForLoad : No. of GPs along boundary
%
%        Output :
% errL2Velocity : The relative error in the L2-norm for the velocity field
% errL2Pressure : The relative error in the L2-norm for the pressure field
%    minElASize : The minimum element area size in the isogeometric mesh
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
%    3i. Initialize the element area size
%
%   3ii. Compute the determinant of the Jacobian to the transformation from the NURBS space (xi-eta) to the integration domain [-1,1]x[-1,1]
%
%  3iii. Loop over all the Quadrature Points
%  ->
%        3iii.1. Compute the coordinates of the quadrature points at the NURBS domain via mapping
%
%        3iii.2. Compute the basis functions on the Gauss Point
%
%        3iii.3. Compute the Cartesian location of the Gauss Point
%
%        3iii.4. Get the element nodal discrete vector
%
%        3iii.5. Compute the numerical solution vector [ux uy p]' on the Gauss Point
%
%        3iii.6. Compute the exact solution on the Gauss Point
%
%        3iii.7. Compute the Jacobian of the transformation from the physical to the parameter space
%
%        3iii.8. Compute the element area on the Gauss Point and sum up the contribution to it
%
%        3iii.9. Compute the square of the L2-norm of the error and the exact solution at the Gauss Point and sum the contribution      
%  <-
%   3iv. Check if the current minimum element area size is smaller than for the previous element
% <-
% 4. Compute the relative error in the L2-norm
%
%% Function main body

%% 0. Read input

% Compute the number of the knots in xi-,eta-direction
mxi = length(Xi);
meta = length(Eta);

% Number of Control Points in xi-direction
nxi = length(CP(:, 1, 1));

% Initialize minimum element area size
minElementSize = 1e4;
minElASize = minElementSize;

% Initialize auxiliary arrays
velocity_analytical = zeros(2, 1);

% Initialize output
errL2Velocity = 0;
L2VelocityExact = 0;
errL2Pressure = 0;
L2PressureExact = 0;

%% 1. Create the element nodal solution vectors
upEl = zeros(mxi - p - 1, meta - q - 1, 3*(p + 1)*(q + 1));
for j = (q + 1):(meta - q - 1)
    for i = (p + 1):(mxi - p - 1)
        k = 1; 
        for c = j - q - 1:j - 1 
            for b = i - p:i
                upEl(i, j, k) = up(3*(c*nxi + b) - 2);
                upEl(i, j, k + 1) = up(3*(c*nxi + b) - 1);
                upEl(i, j, k + 2) = up(3*(c*nxi + b));
                k = k + 3;
            end
        end
    end
end

%% 2. Get the quadrature for the integration of the error
if strcmp(propInt.type,'default')
    nGPXi = p + 1;
    nGPEta = q + 1;
elseif strcmp(propInt.type,'user')
    nGPXi = propInt.xiNGP;
    nGPEta = propInt.etaNGP;
end
[xiGP, xiGW] = getGaussPointsAndWeightsOverUnitDomain(nGPXi);
[etaGP, etaGW] = getGaussPointsAndWeightsOverUnitDomain(nGPEta);

%% 3. Loop over all the elements (knot spans)
for j = (q + 1):(meta - q - 1)   
    for i = (p+1):(mxi-p-1)
        if (Xi(i + 1) ~= Xi(i) && Eta(j + 1) ~= Eta(j))
            %% 3i. Initialize the element area size
            minElementSize = 0;
            
            %% 3ii. Compute the determinant of the Jacobian to the transformation from the NURBS space (xi-eta) to the integration domain [-1,1]x[-1,1] 
            % 
            %         | xi_i+1 - xi_i                    |
            %         | -------------            0       |
            %         |        2                         |
            %  xi,u = |                                  |
            %         |                  eta_j+1 - eta_j |
            %         |        0         --------------- |
            %         |                          2       |
            detJxiu = (Xi(i + 1) - Xi(i))*(Eta(j + 1) - Eta(j))/4;
            
            %% 3iii. Loop over all the Quadrature Points
            for ieta = 1:length(etaGP)
                for ixi =1:length(xiGP)
                    %% 3iii.1. Compute the coordinates of the quadrature points at the NURBS domain via mapping
                    xi = (Xi(i + 1) + Xi(i) + xiGP(ixi)*(Xi(i + 1) - Xi(i)))/2;
                    eta = (Eta(j + 1) + Eta(j) + etaGP(ieta)*(Eta(j + 1) - Eta(j)) )/2;
                    
                    %% 3iii.2. Compute the basis functions on the Gauss Point
                    dR = computeIGABasisFunctionsAndDerivativesForSurface ...
                        (i, p, xi, Xi, j, q, eta, Eta, CP, isNURBS, 1);
                    
                    %% 3iii.3. Compute the Cartesian location of the Gauss Point
                    PCartesian = computeCartesianCoordinatesOfAPointOnBSplineSurface ...
                        (i, p, xi, Xi, j, q, eta, Eta, CP, dR(:,1));
                    x = PCartesian(1);
                    y = PCartesian(2);
                    
                    %% 3iii.4. Get the element nodal discrete vector
                    upActual(:, 1) = upEl(i, j, :);
                    
                    %% 3iii.5. Compute the numerical solution vector [ux uy p]' on the Gauss Point
                    upVct = computeNodalVectorIncompressibleFlow2D(dR(:,1), p, q, upActual);
                    velocity = upVct(1:2, 1);
                    pressure = upVct(3, 1);
                    
                    %% 3iii.6. Compute the exact solution on the Gauss Point
                    velocity_analytical(1, 1) = x^2 * (1 - x)^2 * (2*y - 6*y^2 + 4*y^3);
                    velocity_analytical(2, 1) = -y^2 * (1 - y)^2 * (2*x - 6*x^2 + 4*x^3);
                    pressure_analytical = x*(1-x);
                    
                    %% 3iii.7. Compute the Jacobian of the transformation from the physical to the parameter space
                    Jxxi = zeros(2, 2);
                    k = 0;
                    for c = 0:q
                        for b = 0:p
                            k = k + 1;
                            Jxxi(1, 1) = Jxxi(1, 1) + CP(i -p +b, j - q + c, 1)*dR(k, 2);
                            Jxxi(1, 2) = Jxxi(1, 2) + CP(i -p +b, j - q + c, 2)*dR(k, 2);
                            Jxxi(2, 1) = Jxxi(2, 1) + CP(i -p +b, j - q + c, 1)*dR(k, 3);
                            Jxxi(2, 2) = Jxxi(2, 2) + CP(i -p +b, j - q + c, 2)*dR(k, 3);
                        end
                    end
                    detJxxi = det(Jxxi);
                    
                    %% 3iii.8. Compute the element area on the Gauss Point and sum up the contribution to it
                    minElementSizeOnGP = abs(detJxxi)*abs(detJxiu)*xiGW(ixi)*etaGW(ieta);
                    minElementSize = minElementSize + minElementSizeOnGP;
                    
                    %% 3iii.9. Compute the square of the L2-norm of the error and the exact solution at the Gauss Point and sum the contribution
                    
                    % For the velocity field
                    errL2Velocity = errL2Velocity + norm(velocity - velocity_analytical)^2*minElementSizeOnGP;
                    L2VelocityExact = L2VelocityExact + norm(velocity)^2*minElementSizeOnGP;
                    
                    % For the pressure field
                    errL2Pressure = errL2Pressure + norm(pressure - pressure_analytical)^2*minElementSizeOnGP;
                    L2PressureExact = L2PressureExact + norm(pressure)^2*minElementSizeOnGP;
                end
            end
            
            %% 3iv. Check if the current minimum element area size is smaller than for the previous element
            if minElementSize <= minElASize
                minElASize = minElementSize;
            end
        end
    end
end

%% 4. Compute the relative error in the L2-norm

% Compute the relative error in the L2-norm for the velocity field
errL2Velocity = sqrt(errL2Velocity)/sqrt(L2VelocityExact);

% Compute the relative error in the L2-norm for the pressure field
errL2Pressure = sqrt(errL2Pressure)/sqrt(L2PressureExact);

end