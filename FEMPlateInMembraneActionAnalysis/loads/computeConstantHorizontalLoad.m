function load = computeConstantHorizontalLoad(x, y, z, t, propNBC)
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
% time t. The load is assumed to be constant and horizontal (x-direction).
%
%       Input :
%       x,y,z : The physical location where the load is applied
%           t : The time instance
%     propNBC : Structure defining properties regarding the Neumann
%               boundary conditions
%
%      Output :
%        load :  The load vector [loadx; loady; loadz]
%
%% Function main body
loadAmplitude = 1e3;
load = zeros(3,1);
load(1,1) = loadAmplitude;
if isfield(propNBC, 'tractionVector')
    load = propNBC.tractionVector;
    if isfield(propNBC, 'endTime')
        if isnumeric(propNBC.endTime)
            if t > propNBC.endTime
                load = zeros(3, 1);
            end
        end
    end
end

end
