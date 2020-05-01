function homDOFs = findDofs2D(homDOFs, xi, eta, dir, CP)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% makes fixed supports and adds these to existing ones
% suitable to fix a corner or an edge
%
%   Input :
% homDOFs : previous set of supports
%  xi,eta : region to be supported (e.g. xi = [0 1], eta = [0 1])
%     dir : direction: 1-x, 2-y
%      CP : Control Point coordinates and weights
%
%  Output :
%      rb : new set of supports 
%
%% Function main body

% counter
r = length(homDOFs) + 1;

% number of control points in u,v-direction
numCPs_xi = length(CP(:, 1, 1));
numCPs_eta = length(CP(1, :, 1));

% Iterate and add new supports preserving the old ones
for j = eta(1)*(numCPs_eta - 1) + 1:eta(2)*(numCPs_eta - 1) + 1
    for i = xi(1)*(numCPs_xi - 1) + 1:xi(2)*(numCPs_xi - 1) + 1
        homDOFs(r) = 2*((j - 1)*numCPs_xi + i-1) + dir;
        
        % Round to nearest integer
        homDOFs(r) = round(homDOFs(r));
        
        % Update counter
        r = r + 1;
    end
end

% sort rb and delete double entries
homDOFs = sort(homDOFs);

% Initialize counter
i = 1;

% Loop over the supports and delete double entries if any
while i < length(homDOFs)
    if homDOFs(i) == homDOFs(i + 1)
        homDOFs(i+1) = [];  
        
        % Decrease counter
        i = i - 1;  
    end
    % Update counter
    i = i + 1;
end

end
