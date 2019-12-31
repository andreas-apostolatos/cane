function distribution = generateRandomLatinHypercubeDistribution(meanValue, standardDeviation, noSamples)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns an array of random numbers generated based on the Latin Hypercube 
% distribution.
%
%             Input :
%         meanValue : Mean value of the distribution
% standardDeviation : Standard deviation of the distribution
%         noSamples : Number of samples
%
%            Output :
%      distribution : Random vector of noSamples distributed according to
%                     the Latin Hypercube
%       
%% Function main body
distribution = lhsnorm(meanValue, standardDeviation^2, noSamples);

end