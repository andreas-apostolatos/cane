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
function CPd = computeDisplacedControlPointsForIGAPlateInMembraneAction(CP,dHat)
%% Function documentation
%
% Returns the displaced Control Points corresponding to the isogeometric
% plate in membrane action.
%
%  Input : 
%     CP : control point coordinates and weights before deformation
%   dHat : the Control Point displacement field
%
% Output :
%    CPd : control point coordinates and weights after deformation
%
% Function layout :
%
% 0. Read input
%
% 1. Add the dipslacements into the array CPd iteratively
%
%% Function main body

%% 0. Read input

% Number of control points in u,v-direction
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));

% Initialize counter
k = 1;

% Initialize array of deformed control point location
CPd = zeros(nxi,neta,length(CP(1,1,:)));

%% 1. Add the dipslacements into the array CPd iteratively
for j = 1:neta
    for i = 1:nxi
        % Add odd DOFs to x-component of the CPs
        CPd(i,j,1) = CP(i,j,1) + dHat(k);
        
        % Add even DOFs to y-component of the CPs
        CPd(i,j,2) = CP(i,j,2) + dHat(k+1);
        
        % Preserve the Control Point Weights and the z-component of the CPs 
        CPd(i,j,3:4) = CP(i,j,3:4);
        
        % Update counter 
        k = k + 2;
    end
end

end