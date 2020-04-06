%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _______________________________________________________               %
%   _______________________________________________________               %
%                                                                         %
%   Technische Universität München                                        %
%   Lehrstuhl für Statik, Prof. Dr.-Ing. Kai-Uwe Bletzinger               %
%   _______________________________________________________               %
%   _______________________________________________________               %
%                                                                         %
%                                                                         %
%   Authors                                                               %
%   _______________________________________________________________       %
%                                                                         %
%   Dipl.-Math. Andreas Apostolatos    (andreas.apostolatos@tum.de)       %
%   Dr.-Ing. Roland Wüchner            (wuechner@tum.de)                  %
%   Prof. Dr.-Ing. Kai-Uwe Bletzinger  (kub@tum.de)                       %
%   _______________________________________________________________       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e = computeErrorStokesFlowInUnitSquareUnderBodyForces(p,q,U,V,CP,up,error)
%% Function documentation
%
% Returns the error corresponding to the Finite Element Method applied to
% the Stokes problem in a unit square domain under certain type of body
% forces.
%
%       Input :
%         p,q : Polynomial degrees
%         U,V : Knot vectors in u,v-direction
%          CP : Control point coordinates and weights
%          up : The nodal solution vector on the Control Points
%       error : Structure desicive on the error
%
%      Output :
%           e : The corresponding error
%
% Function layout :
%
%% Function main body

%% 0. Read input

% Compute the control point coordinates for the deformed configuration
CPd = computeDisplacedControlPointsIGAIncompressibleFlow2D(CP,up,error.component);

% Number of knots in u,v-direction
mu = length(U);
mv = length(V);

% Number of Control Points in u,v-direction
nu = length(CPd(:,1,1));
nv = length(CPd(1,:,1));

% Initialize output
e = 0;
eAnalytical = 0;

%% 2. Create the element nodal solution vectors

% Initialization over the array
uEl = zeros(mu-p-1,mv-q-1,3*(p+1)*(q+1));

% Recursive assignment of the elements in the array
for j = (q+1):(mv-q-1)
    for i = (p+1):(mu-p-1)
        k=1; 
        for c = j-q-1:j-1 
            for b = i-p:i
                uEl(i,j,k)   = up(3*(c*nu + b)-2);
                uEl(i,j,k + 1) = up(3*(c*nu + b)-1);
                uEl(i,j,k + 2) = up(3*(c*nu + b));
                k = k + 3;
            end
        end
    end
end

%% 3. Get the quadrature for the integration of the error

% Get the quadrature rule in u-direction
nGPu = error.int.uNGauss;
[GPu,GWu] = gauss(nGPu);

% Get the quadrature rule in v-direction
nGPv = error.int.vNGauss;
[GPv,GWv] = gauss(nGPv);

%% 4. Loop over all the elements (knot spans)
for j = (q+1):(mv-q-1)   
    for i = (p+1):(mu-p-1)
        % check if we are in a non-zero knot span
        if (U(i+1)~=U(i) && V(j+1)~=V(j))
            %% Compute the determinant of the Jacobian to the transformation from the parameter space to the parent domain
            detJxiu = (U(i+1)-U(i))*(V(j+1)-V(j))/4;
            
            %% 4i. Loop over all the Quadrature Points
            for kv = 1:length(GPv)
                for ku =1:length(GPu)
                    %% 4i.1. Compute the coordinates of the quadrature points at the NURBS domain via mapping
                    u = ( U(i+1)+U(i) + GPu(ku)*(U(i+1)-U(i)) )/2;
                    v = ( V(j+1)+V(j) + GPv(kv)*(V(j+1)-V(j)) )/2;
                    
                    %% 4ii.2. Issue the quadrature weights the compute the quadrature weight needed in the multiplication
                    gw = GWu(ku)*GWv(kv);
                    
                    %% 4iii. Get the element nodal transport vector
                    uActual(:,1) = uEl(i,j,:);
                    
                    %% 4iv. Compute the selected component at the Gauss Point
                    upElement = computeNodalVectorIncompressibleFlow2D(i,p,u,U,j,q,v,V,CP,uActual);
                    if error.component==1||error.component==2||error.component==3
                        upFEM = upElement(error.component);
                    elseif error.component==4
                        upFEM = sqrt(upElement(1)^2 + upElement(2)^2);
                    end
                    
                    %% 4v. Compute the Cartesian location of the Gauss Point
                    PCartesian = computePointCartesianCoordinatesOnBSplineSurface(p,i,u,U,q,j,v,V,CP);
                    x = PCartesian(1);
                    y = PCartesian(2);
                    
                    %% Compute the Jacobian to the transformation from the physical to the parameter space
                    [~,dR] = computeNurbsBasisFunctionsAndFirstDerivatives2D(i,p,u,U,j,q,v,V,CP);
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
                    
                    %% 4vi. Compute the analytical values at the Gauss Point location
                    if error.component==1||error.component==2||error.component==4
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
                    
                    %% 4vii. Compute the error at the Gauss Point
                    if strcmp(error.type,'L2')
                        eGP = norm(upFEM - upAnalytical)^2 * gw;
                        eAnalyticalGP = norm(upAnalytical)^2 * gw;
                    end
                    
                    %% 4viii. Add the contribution from the Gauss Point
                    e = e + eGP;
                    eAnalytical = eAnalytical + eAnalyticalGP;
                end
            end
        end
    end
end

% Compute the relative error
e = e/eAnalytical;

end

