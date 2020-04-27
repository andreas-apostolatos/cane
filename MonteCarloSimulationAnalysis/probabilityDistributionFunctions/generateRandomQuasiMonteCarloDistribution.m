function distribution = generateRandomQuasiMonteCarloDistribution ...
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
% Returns an array of random numbers generated based on the Quasi
% Monte-Carlo distribution.
%
%             Input :
%         meanValue : Mean value of the distribution
% standardDeviation : Standard deviation of the distribution
%        numSamples : Number of samples
%
%            Output :
%      distribution : Random vector of noSamples distributed according to
%                     the Quasi Monte-Carlo Simulation
%       
%% Function main body

% Generate necessary variables
p = haltonset(1,'Skip',1e3,'Leap',1e2);
p = scramble(p,'RR2');
haltonvector = net(p,numSamples); % Halton Sequence sampling

% Invert uniform Halton sequence to normal distribution
distribution = norminv(haltonvector, meanValue, standardDeviation);

end
