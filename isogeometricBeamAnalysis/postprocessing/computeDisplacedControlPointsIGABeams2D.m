function CPd = computeDisplacedControlPointsIGABeams2D ...
    (CP, dhat, analysis)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Computes the coordinates of the control points after the deformation for
% both the Bernoulli and the Timoshenko beam settings
% 
%    Input :
%       CP : The control points of the initial configuration
%     dhat : The control variables of the displacement field
% analysis : .type : Bernoulli or Timoshenko beam theory
%
%   Output : 
%
%      CPd : The displaced control points
%
%% Function main body

% Initialize the deformed Control Point locations
CPd = zeros(size(CP,1),4) ;

if strcmp(analysis.type,'Bernoulli')
    for i = 1 : size(CP,1)
        CPd(i,1) = CP(i,1) + dhat(2*i-1);
        CPd(i,2) = CP(i,2) + dhat(2*i);
        CPd(i,3) = CP(i,3);
        CPd(i,4) = CP(i,4);
    end
elseif strcmp(analysis.type,'Timoshenko')
    for i = 1 : size(CP,1)
        CPd(i,1) = CP(i,1) + dhat(3*i-2);
        CPd(i,2) = CP(i,2) + dhat(3*i-1);
        CPd(i,3) = CP(i,3);
        CPd(i,4) = CP(i,4);
    end
end

end
