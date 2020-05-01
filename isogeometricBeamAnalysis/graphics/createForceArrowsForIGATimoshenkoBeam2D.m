function [xF, yF, zF] = createForceArrowsForIGATimoshenkoBeam2D ...
    (CP, F)
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
% to the forces arrows for the case of the isogeometric Timoshenko beam in
% 2D
%
%  Input :
%     CP : Control Point locations
%      F : Complete force vector
%  
% Output :
%     xf : x-coordinates of the points
%     yf : y-coordinates of the points
%     zf : z-coordinates of the points
%
%% Function main body

% Find the number of the force locations 
nF = sum(F~=0);

% Initialize output arrays
xF = zeros(nF,1);
yF = zeros(nF,1);
zF = zeros(nF,1);

% Initialize counters
k = 1;
l = 1;
factor = 1;

% Loop over all the Control Points
for i = 1:length(CP(:,1,1))
    if F(k)~=0
        xF(l,1) = CP(i,1);
        xF(l,2) = factor*(CP(i,1)-F(k))/max(abs(F));
        yF(l,1) = CP(i,2);
        yF(l,2) = CP(i,2);
        zF(l,1) = CP(i,3);
        zF(l,2) = CP(i,3);
        % Update internal counter
        l=l+1;
    end
    if F(k+1)~=0
        xF(l,1) = CP(i,1);
        xF(l,2) = factor*CP(i,1);
        yF(l,1) = CP(i,2);
        yF(l,2) = (CP(i,2)-F(k+1))/max(abs(F));
        zF(l,1) = CP(i,3);
        zF(l,2) = CP(i,3);
        % Update internal counter
        l=l+1;
    end
    if (F(k+2)~=0)
        % This would correspond to directly applying a moment, but if this
        % functions is needed to be extended there are more things to be
        % considered
        
%         xf(l,1) = CP(i,1);
%         xf(l,2) = CP(i,1);
%         yf(l,1) = CP(i,2);
%         yf(l,2) = CP(i,2);
%         zf(l,1) = CP(i,3);
%         zf(l,2) = scal*(CP(i,3)-f(k+2)/max(abs(f)));
% 
%         % Update counter
%         l=l+1;
    end
    
    % Update external counter
    k = k + 3;
end

end
