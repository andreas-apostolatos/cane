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
function T2LC = computeT2LocalCartesianBasis(eA,eLC)
%% Function documentation
%
% Returns the tranformation matrix T for the transformation of a symmetric
% tensor sigma in a Voigt form from the given basis eA to the local 
% Cartesian basis eLC which is defined over a B-Spline surface using the 
% rule eLC1 = G_1/norm(G_1) and eLC2 = G^2/norm(G^2), namely,
%
%   sigmaLC_Voigt = T*sigmaA_Voigt
%
%       Input :
%          eA : The basis where the second order tensor is defined onto
%               given as eA = [eA1 eA2]
%         eLC : The local Cartesian basis on the B-Spline surface
%
%      Output :
%        T2LC : The transformation matrix for a second order symmetric
%               tensor in Voigt form from basis eA to the local Cartesian 
%               basis eLC
%
% Function layout :
%
% 0. Read input
%
% 1. Compute the transformation matrix from from basis eA to the local Cartesian basis eLC
%
%% Function main body

%% 0. Read input

% Initialize output array
T2LC = zeros(3,3);

%% 1. Compute the transformation matrix from from basis eA to the local Cartesian basis eLC
T2LC(1,1) = (eA(:,1)'*eLC(:,1))^2;
T2LC(1,2) = (eA(:,2)'*eLC(:,1))^2;
T2LC(1,3) = 2*(eA(:,1)'*eLC(:,1))*(eA(:,2)'*eLC(:,1));
T2LC(2,1) = (eA(:,1)'*eLC(:,2))^2;
T2LC(2,2) = (eA(:,2)'*eLC(:,2))^2;
T2LC(2,3) = 2*(eA(:,1)'*eLC(:,2))*(eA(:,2)'*eLC(:,2));
T2LC(3,1) = (eA(:,1)'*eLC(:,1))*(eA(:,1)'*eLC(:,2));
T2LC(3,2) = (eA(:,2)'*eLC(:,1))*(eA(:,2)'*eLC(:,2));
T2LC(3,3) = (eA(:,1)'*eLC(:,1))*(eA(:,2)'*eLC(:,2)) + ...
    (eA(:,2)'*eLC(:,1))*(eA(:,1)'*eLC(:,2));

end