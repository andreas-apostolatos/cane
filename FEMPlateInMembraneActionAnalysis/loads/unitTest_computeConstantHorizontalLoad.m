function load = unitTest_computeConstantHorizontalLoad(x,y,z,t)
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
% time t. The load is assumed to be constant and horizontal (x-direction)
% and this function for the computation of the load vector is used in the 
% unit test cases so it is not meant to be modified.
%
%       Input :
%       x,y,z : The physical location where the load is applied
%           t : The time instance
%
%      Output :
%        load :  The load vector [loadx; loady; loadz]
%
%% Function main body

loadAmplitude = -1e0;
load = zeros(3,1);
load(1,1) = loadAmplitude;

end
