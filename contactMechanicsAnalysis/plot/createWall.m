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
function points = createWall(wall)
%% Function documentation
%
% returns coordinates for the contact nodes and the rigid wall marking lines
%
%   Input : 
%    wall : Matrix with the coordinates of two wall determining points
%
%  Output :
%      
%
%% Function main body

% Initialize the output arrays for the rigid wall marking lines
x0(1,2)=0;
x1(1,2)=0;
scaling = 0.05;

% Check whether a wall was defined
if isempty(wall)==0
    
    % change of distance in x and y direction
    DX = wall(2,1)-wall(1,1);
    DY = wall(2,2)-wall(1,2);
    % normalized segment normal
    normal = [-DY,DX];
    normal = normal/norm(normal);
    % Coordinates for the rigid wall marker lines
    x0 = wall(1,:)- normal*scaling;
    x1 = wall(2,:)- normal*scaling; 
end
% output result
points = [x0;x1];

end