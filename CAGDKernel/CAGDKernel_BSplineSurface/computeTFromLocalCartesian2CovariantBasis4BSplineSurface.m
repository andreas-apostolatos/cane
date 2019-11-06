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
function TLC2Cov = computeTFromLocalCartesian2CovariantBasis4BSplineSurface(eLC,GContravariant)
%% Function documentation
%
% Returns the tranformation matrix TLC2Cov which transform a second order
% tensor in Voigt notation from the local Cartesian basis to the covariant
% basis namely:
%
% |n^11|           |n^LC_11|  
% |n^22| = TLC2Cov |n^LC_22|
% |n^12|           |n^LC_12|
%
%          Input :
%            eLC : = [eHat1 eHat2] the local Cartesian basis
% GContravariant : = [GContra1 GContra2 GContra3] the contravariant base 
%                  vectors of the configuration
%
%         Output :
%        TLC2Cov : The transformation matrix from the contravariant 
%                  to the local Cartesian basis
%
% Function layout :
%
% 0. Read input
%
% 1. Compute the transformation matrix from the contravariant to local cartesian basis
%   
%% Function main body

%% 0. Read input

% Initialize output array
TLC2Cov = zeros(3,3);

doubleFactor = 1;

% Compute the entries of the 
TLC2Cov(1,1) = (GContravariant(:,1)'*eLC(:,1))^2;
TLC2Cov(1,2) = (GContravariant(:,1)'*eLC(:,2))^2;
TLC2Cov(1,3) = 2*(GContravariant(:,1)'*eLC(:,1))*(GContravariant(:,1)'*eLC(:,2));
TLC2Cov(2,1) = (GContravariant(:,2)'*eLC(:,1))^2;
TLC2Cov(2,2) = (GContravariant(:,2)'*eLC(:,2))^2;
TLC2Cov(2,3) = 2*(GContravariant(:,2)'*eLC(:,2))*(GContravariant(:,2)'*eLC(:,1));
TLC2Cov(3,1) = doubleFactor*((GContravariant(:,1)'*eLC(:,1))*(GContravariant(:,2)'*eLC(:,1)));
TLC2Cov(3,2) = doubleFactor*((GContravariant(:,1)'*eLC(:,2))*(GContravariant(:,2)'*eLC(:,2)));
TLC2Cov(3,3) = doubleFactor*((GContravariant(:,1)'*eLC(:,1))*(GContravariant(:,2)'*eLC(:,2)) + ...
    (GContravariant(:,1)'*eLC(:,2))*(GContravariant(:,2)'*eLC(:,1)));

end

