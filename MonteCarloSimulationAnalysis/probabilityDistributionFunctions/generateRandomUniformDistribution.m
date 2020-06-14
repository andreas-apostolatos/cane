function distribution = generateRandomUniformDistribution ...
    (meanValue, standardDeviation, numSamples)
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
%        numSamples : Number of samples
%
%            Output :
%      distribution : Random vector of noSamples uniformly distributed
%       
%% Function main body
distribution = 2*standardDeviation*rand(numSamples,1) + meanValue - standardDeviation;

end
