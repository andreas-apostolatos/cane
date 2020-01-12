%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Lehrstuhl für Statik, Prof. Dr.-Ing. K.-U. Bletzinger                 %
%   _____________________________________________________                 %
%                                                                         %
%   Authors                                                                %
%   _______                                                               %
%                                                                         %
%   Dr.-Ing. Roland Wüchner                                               %
%   Dipl.-Math. Andreas Apostolatos (andreas.apostolatos@tum.de)          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sol= createWall(wall)
%% Function documentation
%
% returns coordinates for the contact nodes and the rigid wall marking lines
%
%   Input : 
%    wall : Matrix with the coordinates of two wall determining points
%
%  Output :
%      Structure sol with the following fields:
%
%      xw : x-coordinates of 30 equidistant points on the wall (starting position for 
%           rigid wall marker lines)
%      yw : y-coordinates of 30 equidistant points on the wall (starting position for 
%           rigid wall marker lines)
%
%
% Function layout :
%
%  2.Create rigid wall
%
%  3.Output results
%
%% Function main body

%% 2.Create rigid wall
% Assign a value as number of marking lines on plotting the rigid wall 
number_marking_lines = 30;

% Initialize the output arrays for the rigid wall marking lines
xw(number_marking_lines,2)=0;
yw(number_marking_lines,2)=0;

% Check whether a wall was defined
if isempty(wall)==0

% Direction vector of the wall line
dir=[wall(2,1)-wall(1,1);wall(2,2)-wall(1,2)];

%rotation Matrix for rigid wall marker lines
a=-3*pi/4;
rotation=[cos(a),-sin(a);sin(a),cos(a)];
rot_dir=rotation*dir;
rot_dir=rot_dir./20;

%Coordinates for the  Matrix for rigid wall marker lines
for j=1:number_marking_lines

        xw(j,1)=wall(1,1)+j*1/number_marking_lines*dir(1);
        xw(j,2)= xw(j,1)+rot_dir(1);

        yw(j,1)=wall(1,2)+j*1/number_marking_lines*dir(2);
        yw(j,2)= yw(j,1)+rot_dir(2);
    
        
end
end
%% 3.Output results
sol=struct('xw',xw,'yw',yw);

end