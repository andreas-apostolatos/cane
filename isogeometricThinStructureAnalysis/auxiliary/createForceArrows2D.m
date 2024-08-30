function [xf, yf, zf] = createForceArrows2D(CP, f)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the coordinates of the points for the creation of cooresponding 
% to the forces arrows
%
%    Input :
%       CP : Control Point coordinates
%        f : Complete force vector
%
%   Output :
% xf,yf,zf : x-,y- and z- coordinates of the points where forces are
%            applied
%
%% Function main body

% Find the number of the force locations 
numF = sum(f ~= 0);

% Initialize output arrays
xf = zeros(numF, 2);
yf = zeros(numF, 2);
zf = zeros(numF, 2);

% Initialize counters
counterX = 1;
counterY = 1;

% Loop over all the Control Points
for j = 1:length(CP(1, :, 1))
    for i = 1:length(CP(:, 1, 1))
        if f(counterX) ~= 0
            xf(counterY, 1) = CP(i, j, 1);
            xf(counterY, 2) = CP(i, j, 1) - f(counterX)/max(abs(f));
            yf(counterY, 1) = CP(i, j, 2);
            yf(counterY, 2) = CP(i, j, 2);
            zf(counterY, 1) = CP(i, j, 3);
            zf(counterY, 2) = CP(i, j, 3);
            
            % Update internal counter
            counterY = counterY + 1;
        end
        if f(counterX + 1) ~= 0
            xf(counterY, 1) = CP(i, j, 1);
            xf(counterY, 2) = CP(i, j, 1);
            yf(counterY, 1) = CP(i, j, 2);
            yf(counterY, 2) = CP(i, j, 2) - f(counterX + 1)/max(abs(f));
            zf(counterY, 1) = CP(i, j, 3);
            zf(counterY, 2) = CP(i, j, 3);
            
            % Update internal counter
            counterY = counterY + 1;
        end
        
        % Update external counter
        counterX = counterX + 2;
    end
end

end