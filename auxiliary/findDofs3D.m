function homDOFs = findDofs3D(homDOFs,xi,eta,dir,CP)
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
%     dir : direction: 1-x, 2-y, 3-z
%
%  Output :
%      rb : new set of supports 
%
%% Function main body

% initialize counter
r = length(homDOFs) + 1;

% number of control points in xi,eta-direction
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));

% Iterate and add new supports preserving the old ones
for j = eta(1)*(neta-1)+1:eta(2)*(neta-1)+1
    for i = xi(1)*(nxi-1)+1:xi(2)*(nxi-1)+1
        % Find the corresponding DOF
        homDOFs(r) = 3*((j-1)*nxi + i-1) + dir;
        
        % Update counter
        r = r + 1;
    end
end

% sort homDOFs and delete double entries
homDOFs = sort(homDOFs);
homDOFs = unique(homDOFs);

end
