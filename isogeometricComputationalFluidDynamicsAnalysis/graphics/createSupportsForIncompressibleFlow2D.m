function [xs, ys, zs] = createSupportsForIncompressibleFlow2D ...
    (CP, homDOFs)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns coordinates for the triangles at the support locations for a 2D
% incompressible flow using isogeometric analysis.
%
%   Input : 
%      CP : Control Point coordinates
% homDOFs : Array containing information on the supports
%
%  Output :
%      xs : x-coordinates of the support triangle vertices
%      ys : y-coordinates of the support triangle vertices
%      zs : z-coordinates of the support triangle vertices
%
%% Function main body

% Total number of Control Points
nu = length(CP(:,1,1));

% scaling factors for the support triangles
up = max(max(max(max(CP))));
lo = min(min(min(min(CP))));

% Average the factor with respect to the maximum and minimum values 
fac = (up-lo)/5;

% Sampling points for the rotational supports
numEval = 3;

% Get a step size in the 0 to 2pi interval
noOfCircularSegments = 4;
step = 2*pi/noOfCircularSegments/numEval;

% Initialize a counter for the supports
counterSupports = 1;

% Initialize the support Cartesian coordinates
xs = [];
ys = [];
zs = [];

for l = 1:length(homDOFs)
    % get the corresponding Control Point number p and indices CP(i,j)
    p=ceil(homDOFs(l)/3);
    j=ceil(p/nu);
    i=p-(j-1)*nu;
    dir=homDOFs(l)-((j-1)*nu+i-1)*3;
    
    if (dir==1)              % (x-support)
        xs(counterSupports,1)=CP(i,j,1);
        xs(counterSupports,2)=CP(i,j,1) - 0.1732*fac;
        xs(counterSupports,3)=CP(i,j,1) - 0.1732*fac;
        xs(counterSupports,4)=xs(counterSupports,1);
        ys(counterSupports,1)=CP(i,j,2);
        ys(counterSupports,2)=CP(i,j,2) + 0.1*fac;
        ys(counterSupports,3)=CP(i,j,2) - 0.1*fac;
        ys(counterSupports,4)=ys(counterSupports,1);
        zs(counterSupports,1:4)=CP(i,j,3);
        
        % Update the support counter
        counterSupports = counterSupports + 1;
    elseif (dir==2)          % (y-support)
        xs(counterSupports,1)=CP(i,j,1);
        xs(counterSupports,2)=CP(i,j,1) - 0.1*fac;
        xs(counterSupports,3)=CP(i,j,1) + 0.1*fac;
        xs(counterSupports,4)=xs(counterSupports,1);
        ys(counterSupports,1)=CP(i,j,2);
        ys(counterSupports,2)=CP(i,j,2) - 0.1732*fac;
        ys(counterSupports,3)=CP(i,j,2) - 0.1732*fac;
        ys(counterSupports,4)=ys(counterSupports,1);
        zs(counterSupports,1:4)=CP(i,j,3);
        
        % Update the support counter
        counterSupports = counterSupports + 1;
    elseif (dir==3)          % (rotational-support)
        % Initialize the coordinate on the circle
        uOnCircle = 0;
        
        % Create the circular segments
        for g=1:noOfCircularSegments
            for r=1:numEval+1
                xs(counterSupports,r) = CP(i,j,1) + 0.1232*fac*sin(uOnCircle);
                ys(counterSupports,r) = CP(i,j,2) + 0.1232*fac*cos(uOnCircle);
                zs(counterSupports,r) = CP(i,j,3);
            
                % Update coordinate on circle
                if r~=numEval+1
                    uOnCircle = uOnCircle + step;
                end
            end
            
             % Update the support counter
            counterSupports = counterSupports + noOfCircularSegments;
        end
    end
end

end
