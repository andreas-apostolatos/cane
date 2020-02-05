function distribution = generateRandomUniformDistribution(meanValue, standardDeviation, noSamples)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns an array of random numbers generated based on the uniform
% distribution.
%
%             Input :
%         meanValue : Mean value of the distribution
% standardDeviation : Standard deviation of the distribution
%         noSamples : Number of samples
%
%            Output :
%      distribution : Random vector of noSamples uniformly distributed
%       
%% Function main body
intervalSize = 0.5*standardDeviation*sqrt(12);
distribution = 2*intervalSize*rand(noSamples,1) + meanValue - intervalSize;

end