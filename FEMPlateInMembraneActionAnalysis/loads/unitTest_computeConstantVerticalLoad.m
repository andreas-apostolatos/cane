function load = unitTest_computeConstantVerticalLoad(x, y, z, t, propNBC)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the applied load vector at the physical location x,y,z and at
% time t. The load is assumed to be constant and vertical (y-direction).
%
%       Input :
%       x,y,z : The physical location where the load is applied
%           t : The time instance
%     propNBC : Structure containing information on the Neumann boundary
%               conditions
%
%      Output :
%        load :  The load vector [loadx; loady; loadz]
%
%% Function main body
loadAmplitude = -1e2;
load = zeros(3,1);
load(2,1) = loadAmplitude;

end
