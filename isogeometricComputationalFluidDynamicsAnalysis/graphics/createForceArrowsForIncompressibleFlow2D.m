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
function [xf,yf,zf] = createForceArrowsForIncompressibleFlow2D(CP,f)
%% Function documentation
%
% Returns the coordinates of the points for the creation of cooresponding 
% to the forces arrows (flux arrows) for the case of a 2D incompressible
% flow problem.
%
%  Input :
%     CP : Control Point locations
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
k=1;
l=1;

% Loop over all the Control Points
for j = 1:length(CP(1,:,1))
    for i = 1:length(CP(:,1,1))
        if f(k)~=0
            xf(l,1) = CP(i,j,1);
            xf(l,2) = CP(i,j,1)-f(k)/max(abs(f));
            yf(l,1) = CP(i,j,2);
            yf(l,2) = CP(i,j,2);
            zf(l,1) = CP(i,j,3);
            zf(l,2) = CP(i,j,3);
            % Update internal counter
            l=l+1;
        end
        if f(k+1)~=0
            xf(l,1) = CP(i,j,1);
            xf(l,2) = CP(i,j,1);
            yf(l,1) = CP(i,j,2);
            yf(l,2) = CP(i,j,2)-f(k+1)/max(abs(f));
            zf(l,1) = CP(i,j,3);
            zf(l,2) = CP(i,j,3);
            % Update internal counter
            l=l+1;
        end
        % Update external counter
        k=k+3;
    end
end

end