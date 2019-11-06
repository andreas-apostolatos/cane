function [dExact,betaExact] = computeExactDispl4BernoulliCantileverBeamInUniformPressure(P,problemSettings)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the exact displacement field at the physical location P of a
% cantilever Euler-Bernoulli beam which is subject to a uniform pressure
% load and which's two ends are in the physical locations (0,0) - (0,L).
%
%           Input :
%               P : The physical location where to compute the exact 
%                   solution to the problem
% problemSettings : Parameters which are related to the benchmark
%                   problem, such as dimensions, load amplitude etc.
%
%          Output :
%          dExact : The exact displacement solution of the benchmark 
%                   problem at P
%       betaExact : The exact rotation solution of the benchmark problem
%                   at P (This is a dummy output for this function)
%
% Function layout :
%
% 0. Read input
%
% 1. Compute the exact solution at the given Cartesian location
%
%% Function main body

%% 0. Read input

% Re-assign problem settings
L = problemSettings.Length;
pLoad = problemSettings.pressure;
EYoung = problemSettings.EYoung;
I = problemSettings.I;
x = P(1,1);

% Initialize output array
dExact = zeros(2,1);

%% 1. Compute the exact solution at the given Cartesian location
dExact(1,1) = 0;
dExact(2,1) = pLoad*x^2*(6*L^2 - 4*L*x + x^2)/24/EYoung/I;
betaExact = 'undefined';

end

