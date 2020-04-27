function [dA3KommaAlpha, dAKommaAlpha] = ...
    computeParametricDrvsSurfaceNormalOnBSplineSurface ...
    (GCovariant, dGCovariant, A3, dA)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the parametric derivatives along xi- and eta- directions of the 
% surface normal of a B-Spline surface given the base vectors at the
% parametric location of interest.
%
%         Input :
%    GCovariant : The covariant base vectors of the curvilinear 
%                 coordinate system
%   dGCovariant : The derivatives of the covariant base vectors
%            A3 : The unit normal to the NURBS surface base vector
%            dA : The differential element area
%
%        Output :
% dA3KommaAlpha : The parametric derivatives of the surface normal vector 
%                 = [dA3/dxi dA3/deta]
%  dAKommaAlpha : The parametric derivatives of the element surface area
%
% Function layout :
%
% 0. Read input
%
% 1. Compute the not normalized normal to the surface vector
%
% 2. Compute parametric derivatives of the not normalized surface normal vector
%
% 3. Compute the derivatives of the normalized surface normal vector as dA3 = [dA3/d(theta^1) dA3/d(theta^1)]
%
%% Function main body

%% 0. Read input

% Initialize output
dA3KommaAlpha = zeros(3, 2);
dAKommaAlpha = zeros(2, 1);

%% 1. Compute the not normalized normal to the surface vector
A3Tilde = A3*dA;

%% 2. Compute parametric derivatives of the not normalized surface normal vector
A1TimesA2KommaAlpha = zeros(3, 2);
A1TimesA2KommaAlpha(:, 1) = ...
    cross(dGCovariant(:, 1), GCovariant(:, 2)) + cross(GCovariant(:, 1), dGCovariant(:, 3));
A1TimesA2KommaAlpha(:, 2) = ...
    cross(dGCovariant(:, 3), GCovariant(:, 2)) + cross(GCovariant(:, 1), dGCovariant(:, 2));

%% 3. Compute the derivatives of the normalized surface normal vector as dA3 = [dA3/d(theta^1) dA3/d(theta^1)]
for alpha = 1:2
    % Parametric derivatives of the element surface area
    dAKommaAlpha(alpha,1) = (1/dA)*A3Tilde'*A1TimesA2KommaAlpha(:,alpha);
    
    % Parametric derivatives of the normal to the surface base vector
    dA3KommaAlpha(:,alpha) = (1/dA)*A1TimesA2KommaAlpha(:,alpha) - ...
        (1/dA^3)*(A3Tilde'*A1TimesA2KommaAlpha(:,alpha))*A3Tilde;
end

end
