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
function TLC2Contra = computeTFromLocalCartesian2ContravariantBasis4BSplineSurface(eLC,GCovariant)
%% Function documentation
%
% Returns the tranformation matrix TLC2Cov which transform a second order
% tensor in Voigt notation from the local Cartesian basis to the
% contravariant basis namely:
%
% |n_11|              |n^LC_11|  
% |n_22| = TLC2Contra |n^LC_22|
% |n_12|              |n^LC_12|
%
%      Input :
%        eLC : = [eHat1 eHat2] the local Cartesian basis
% GCovariant : = [GCov1 GCov2 GCov3] the contravariant base 
%              vectors of the configuration
%
%     Output :
% TLC2Contra : The transformation matrix from the contravariant 
%              to the local Cartesian basis
%   
%% Function main body

% Initialize output array
TLC2Contra = zeros(3,3);

% Compute the entries of the 
TLC2Contra(1,1) = (GCovariant(:,1)'*eLC(:,1))^2;
TLC2Contra(1,2) = (GCovariant(:,1)'*eLC(:,2))^2;
TLC2Contra(1,3) = 2*(GCovariant(:,1)'*eLC(:,1))*(GCovariant(:,1)'*eLC(:,2));
TLC2Contra(2,1) = (GCovariant(:,2)'*eLC(:,1))^2;
TLC2Contra(2,2) = (GCovariant(:,2)'*eLC(:,2))^2;
TLC2Contra(2,3) = 2*(GCovariant(:,2)'*eLC(:,2))*(GCovariant(:,2)'*eLC(:,1));
TLC2Contra(3,1) = ((GCovariant(:,1)'*eLC(:,1))*(GCovariant(:,2)'*eLC(:,1)));
TLC2Contra(3,2) = ((GCovariant(:,1)'*eLC(:,2))*(GCovariant(:,2)'*eLC(:,2)));
TLC2Contra(3,3) = ((GCovariant(:,1)'*eLC(:,1))*(GCovariant(:,2)'*eLC(:,2)) + ...
    (GCovariant(:,1)'*eLC(:,2))*(GCovariant(:,2)'*eLC(:,1)));

end

