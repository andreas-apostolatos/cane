function rb = findDofs2D(rb,xi,eta,dir,CP)
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
%  Input :
%     rb : previous set of supports
% xi,eta : region to be supported (e.g. xi = [0 1], eta = [0 1])
%    dir : direction: 1-x, 2-y
%     CP : Control Point coordinates and weights
%
% Output :
%     rb : new set of supports 
%
%% Function main body

% counter
r = length(rb) + 1;

% number of control points in u,v-direction
nu = length(CP(:,1,1));
nv = length(CP(1,:,1));

% Iterate and add new supports preserving the old ones
for j = eta(1)*(nv-1)+1:eta(2)*(nv-1)+1
    for i = xi(1)*(nu-1)+1:xi(2)*(nu-1)+1
        rb(r) = 2*((j-1)*nu + i-1) + dir;
        
        % Round to nearest integer
        rb(r) = round(rb(r));
        
        % Update counter
        r = r + 1;
    end
end

% sort rb and delete double entries
rb = sort(rb);

% Initialize counter
i=1;

% Loop over the supports and delete double entries if any
while i < length(rb)
    if rb(i) == rb(i+1)
        rb(i+1) = [];  
        
        % Decrease counter
        i = i - 1;  
    end
    % Update counter
    i = i + 1;
end

end
