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
function eLC = computeLocalCartesianBasis4BSplineSurface(GCovariant,GContravariant)
%% Function documentation
%
% Returns the local Cartesian basis [eHat1 eHa2] over a B-Spline surface.
% This function is used in the computation of the stiffness to the
% isogeometric Kirchhoff-Love shell
%
%           Input :
%      GCovariant : = [GCov1 GCov2 GCov3] the covariant base vectors of the
%                   configuration
%  GContravariant : = [GContra1 GContra2 GContra3] the covariant base 
%                   vectors of the configuration
%
%          Output :
%             eLC : The local Cartesian base vectors
%
% Function layout :
%
% 0. Read input
%
% 1. Compute the local Cartesian basis
%
%% Function main body

%% 0. Read input

% Initialize output array
eLC = zeros(3,2);

%% 1. Compute the local Cartesian basis

% Normalize the contravariant base vectors (G_1_hat and G^2_hat will be 
% used as the local Cartesian basis since they are orthogonal to each other
% yet both being tangent to the surface)

% As first use the covariant base vector GCov1
normGCov1 = norm(GCovariant(:,1));
eLC(:,1) = GCovariant(:,1)/normGCov1;

% As second use the contravariant base vector GContra2
normGContra2 = norm(GContravariant(:,2));
eLC(:,2) = GContravariant(:,2)/normGContra2;

% As third base vector use the normalized normal to the surface vector G3
% eLC(:,3) = GCovariant(:,3);

end

