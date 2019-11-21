function [xf,yf,zf] = createForceArrows(nodes,f)
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
%  Input :
%  nodes : The coordinates of the nodes
%      f : Complete force vector
%  
% Output :
%     xf : x-coordinates of the points
%     yf : y-coordinates of the points
%     zf : z-coordinates of the points
%
%% Function main body

% Find the number of the force locations 
nf = sum(f~=0);

% Initialize output arrays
xf = zeros(nf,2);
yf = zeros(nf,2);
zf = zeros(nf,2);

% Initialize counters
l=1;
k=1;

% Loop over all the Control Points
for i=1:length(nodes(:,1))
    
    % Get the corresponding Control Point number p and indices node(k)
%     p=ceil(i/2);
    
    if f(k)~=0
        xf(l,1) = nodes(i,1);
        xf(l,2) = nodes(i,1)-f(k)/max(abs(f));
        yf(l,1) = nodes(i,2);
        yf(l,2) = nodes(i,2);
        zf(l,1) = nodes(i,3);
        zf(l,2) = nodes(i,3);
        
        % Update internal counter
        l=l+1;
    end
    if f(k+1)~=0
        xf(l,1) = nodes(i,1);
        xf(l,2) = nodes(i,1);
        yf(l,1) = nodes(i,2);
        yf(l,2) = nodes(i,2)-f(k+1)/max(abs(f));
        zf(l,1) = nodes(i,3);
        zf(l,2) = nodes(i,3);
        
        % Update internal counter
        l=l+1;
    end
    
    % Update external counter
    k=k+2;
end

end
