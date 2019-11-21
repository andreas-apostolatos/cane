function bF = computeConstantVerticalBodyForceVct(x,y,z,t)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the load corresponding to a constant body force at the given
% location x,y,z and at time t.
%
%       Input :
%       x,y,z : The physical location x,y,z
%           t : the time instance
%
%      Output :
%          bF : The body force vector bF = [0; 0; 0]
%
%% Function main body

bF = zeros(length(x),1,3);
% bF = zeros(3,1);

end

