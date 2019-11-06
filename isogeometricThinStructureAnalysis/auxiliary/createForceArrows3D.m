function [xf,yf,zf] = createForceArrows3D(CP,Fl)
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
%     CP : Control Point locations
%     Fl : Complete force vector
%  
% Output :
%     xf : x-coordinates of the points
%     yf : y-coordinates of the points
%     zf : z-coordinates of the points
%
%% Function main body

% Find the number of the force locations 
nf = sum(Fl~=0);

% Initialize the vectors
xf=zeros(nf,1);
yf=zeros(nf,1);
zf=zeros(nf,1);

% Initialize counters
k = 1;
l = 1;
scal = 1;

% Loop over the Control Points
for j = 1:length(CP(1,:,1))
    for i = 1:length(CP(:,1,1))
        if (Fl(k)~=0)
            xf(l,1) = CP(i,j,1);
            xf(l,2) = scal*(CP(i,j,1)-Fl(k)/max(abs(Fl)));
            yf(l,1) = CP(i,j,2);
            yf(l,2) = CP(i,j,2);
            zf(l,1) = CP(i,j,3);
            zf(l,2) = CP(i,j,3);
            
            % Update counter
            l=l+1;
        end
        if (Fl(k+1)~=0)
            xf(l,1) = CP(i,j,1);
            xf(l,2) = CP(i,j,1);
            yf(l,1) = CP(i,j,2);
            yf(l,2) = scal*(CP(i,j,2)-Fl(k+1)/max(abs(Fl)));
            zf(l,1) = CP(i,j,3);
            zf(l,2) = CP(i,j,3);
            
            % Update counter
            l=l+1;
        end
        if (Fl(k+2)~=0)
            xf(l,1) = CP(i,j,1);
            xf(l,2) = CP(i,j,1);
            yf(l,1) = CP(i,j,2);
            yf(l,2) = CP(i,j,2);
            zf(l,1) = CP(i,j,3);
            zf(l,2) = scal*(CP(i,j,3)-Fl(k+2)/max(abs(Fl)));
            
            % Update counter
            l=l+1;
        end
        
        % Update counter
        k=k+3;
    end
end

  
end
