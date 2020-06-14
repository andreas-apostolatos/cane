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
% Returns the applied flux corresponding to the thermal conduction problem
% dependent on the geometric location x, y, t and time t.
%
%       Input :
%       x,y,z : The physical location where the load is applied
%           t : The time instance
%     propNBC : User-defined parameters for the load,
%                  .flux : Externally applied flux
%
%      Output :
%        flux :  The applied flux
%
%% Function main body
flux = propNBC.flux;

end
