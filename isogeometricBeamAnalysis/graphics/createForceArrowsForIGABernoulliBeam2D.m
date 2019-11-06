function [xf,yf,zf] = createForceArrowsForIGABernoulliBeam2D(CP,F)
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
% to the forces arrows for the isogeometric Benroulli beam in 2D
%
%  Input :
%     CP : Control Point locations
%      F : The force vector
%  
% Output :
%     xf : x-coordinates of the points
%     yf : y-coordinates of the points
%     zf : z-coordinates of the points
%
%% Function main body

% Find the number of the force locations 
nf = sum(F~=0);

% Initialize output arrays
xf = zeros(nf,2);
yf = zeros(nf,2);
zf = zeros(nf,2);

% Initialize counters
k = 1;
l = 1;

% Loop over all the Control Points
for i = 1:length(CP(:,1,1))
    if F(k)~=0
        xf(l,1) = CP(i,1);
        xf(l,2) = CP(i,1)-F(k)/max(abs(F));
        yf(l,1) = CP(i,2);
        yf(l,2) = CP(i,2);
        zf(l,1) = CP(i,3);
        zf(l,2) = CP(i,3);
        
        % Update internal counter
        l = l + 1;
    end
    if F(k+1)~=0
        xf(l,1) = CP(i,1);
        xf(l,2) = CP(i,1);
        yf(l,1) = CP(i,2);
        yf(l,2) = CP(i,2)-F(k+1)/max(abs(F));
        zf(l,1) = CP(i,3);
        zf(l,2) = CP(i,3);
        
        % Update internal counter
        l = l + 1;
    end
    
    % Update external counter
    k = k + 2;
end

end
