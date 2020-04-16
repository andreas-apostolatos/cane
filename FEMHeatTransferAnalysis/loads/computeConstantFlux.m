function load = computeConstantLoad(x, y, z, t)
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
%     propNBC : User-defined parameters for the load,
%                  .tractionLoadVct : Externally applied traction vector
%
%      Output :
%        load :  The load vector [loadx; loady; loadz]
%
%% Function main body
load = [1e5; 0; 0];

end

