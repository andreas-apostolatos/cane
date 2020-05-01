function flux = computeConstantFlux(x, y, z, t, propNBC)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%                  Marko Leskovar
%
%% Function documentation
%
% Returns the applied fux at the given chorochronical location.
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
flux = propNBC.flux;

end