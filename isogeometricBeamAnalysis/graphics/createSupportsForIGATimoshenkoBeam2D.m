function [xs,ys,zs] = createSupportsForIGATimoshenkoBeam2D(CP,rb)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns coordinates for the triangles at the support locations which are
% pescribed displacements and the circular arcs to those locations where
% rotations are prescribed for the isogeometric Timoshenko beam in 2D.
%
%   Input : 
%      CP : Control Point coordinates
%      rd : Array containing information on the supports
%
%  Output :
%      xs : x-coordinates of the support triangle vertices
%      ys : y-coordinates of the support triangle vertices
%      zs : z-coordinates of the support triangle vertices
%
%% Function main body

% Total number of Control Points
nu = length(CP(:,1));

% scaling factors for the support triangles
maximum = max(max(max(max(CP))));
minimum = min(min(min(min(CP))));

% Average the factor with respect to the maximum and minimum values 
factor = (maximum-minimum)/5;

% Sampling points for the rotational supports
numEval = 3;

% Get a step size in the 0 to 2pi interval
noOfCircularSegments = 4;
step = 2*pi/noOfCircularSegments/numEval;

% Initialize a counter for the supports
counterSupports = 1;

for l = 1:length(rb)
    % get the corresponding Control Point number p and indices CP(i,j)
    p=ceil(rb(l)/3);
    j=ceil(p/nu);
    i=p-(j-1)*nu;
    dir=rb(l)-((j-1)*nu+i-1)*3;
    
    if (dir==1)              % (x-support)
        xs(counterSupports,1)=CP(i,1);
        xs(counterSupports,2)=CP(i,1) - 0.1732*factor;
        xs(counterSupports,3)=CP(i,1) - 0.1732*factor;
        xs(counterSupports,4)=xs(counterSupports,1);
        ys(counterSupports,1)=CP(i,2);
        ys(counterSupports,2)=CP(i,2) + 0.1*factor;
        ys(counterSupports,3)=CP(i,2) - 0.1*factor;
        ys(counterSupports,4)=ys(counterSupports,1);
        zs(counterSupports,1:4)=CP(i,3);
        
        % Update the support counter
        counterSupports = counterSupports + 1;
    elseif (dir==2)          % (y-support)
        xs(counterSupports,1)=CP(i,1);
        xs(counterSupports,2)=CP(i,1) - 0.1*factor;
        xs(counterSupports,3)=CP(i,1) + 0.1*factor;
        xs(counterSupports,4)=xs(counterSupports,1);
        ys(counterSupports,1)=CP(i,2);
        ys(counterSupports,2)=CP(i,2) - 0.1732*factor;
        ys(counterSupports,3)=CP(i,2) - 0.1732*factor;
        ys(counterSupports,4)=ys(counterSupports,1);
        zs(counterSupports,1:4)=CP(i,3);
        
        % Update the support counter
        counterSupports = counterSupports + 1;
    elseif (dir==3)          % (rotational-support)
        % Initialize the coordinate on the circle
        uOnCircle = 0;
        
        % Create the circular segments
        for g=1:noOfCircularSegments
            for r=1:numEval+1
                xs(counterSupports,r) = CP(i,1) + 0.1232*factor*sin(uOnCircle);
                ys(counterSupports,r) = CP(i,2) + 0.1232*factor*cos(uOnCircle);
                zs(counterSupports,r) = CP(i,3);
            
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
