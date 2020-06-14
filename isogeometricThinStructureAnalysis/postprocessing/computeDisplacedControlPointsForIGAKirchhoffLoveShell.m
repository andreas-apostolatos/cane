function CPd = computeDisplacedControlPointsForIGAKirchhoffLoveShell ...
    (CP, dHat)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the displaced Control Points corresponding to the isogeometric
% Kirchhoff-Love shell.
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

% Number of control points in xi,eta-direction
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));

% Initialize counter
k = 1;

% Initialize array of deformed control point location
CPd = zeros(nxi,neta,length(CP(1,1,:)));

%% 1. Add the dipslacements into the array CPd iteratively
for j = 1:neta
    for i = 1:nxi
        % Add the x-component of the CPs
        CPd(i,j,1) = CP(i,j,1) + dHat(k);
        
        % Add the y-component of the CPs
        CPd(i,j,2) = CP(i,j,2) + dHat(k+1);
        
        % Add the z-component of the CPs
        CPd(i,j,3) = CP(i,j,3) + dHat(k+2);
        
        % Preserve the Control Point Weights of the CPs 
        CPd(i,j,4) = CP(i,j,4);
        
        % Update counter 
        k = k + 3;
    end
end

end
