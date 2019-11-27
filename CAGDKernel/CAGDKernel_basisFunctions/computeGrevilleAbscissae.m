function grevilleAbscissae = computeGrevilleAbscissae(indexCP,p,Xi)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the Greville Abscissae parametric location corresponding to the 
% indexCP-th Control Point, i.e. the location where Control Point with
% index indexCP has close to its maximum influence or where the
% corresponding B-Spline basis function attains its maximum value.
%
%             Input :
%           indexCP : Knot span index
%                 p : Polynomial order of the 1D B-Spline basis
%                Xi : Knot vector of the 1D B-Spline basis
%
%            Output :
% grevilleAbscissae : The Greville Abscissae parametric location
%
% Function layout :
%
% 0. Read input
%
% 1. Compute the Greville Abscissae
%
% 2. Check output
%
%% Function main body

%% 0. Read input

% Initialize output
grevilleAbscissae = 0;

% Define tolerance
eps = 1e-6;

%% 1. Compute the Greville Abscissae
for i = 1:p
    grevilleAbscissae = grevilleAbscissae + Xi(indexCP + i);
end
grevilleAbscissae = grevilleAbscissae/p;

%% 2. Check output
if grevilleAbscissae < Xi(1) - eps || grevilleAbscissae > Xi(end) + eps
    error('Greville Abscissae found outside the parametric space of the curve');
end
if grevilleAbscissae < Xi(1) && grevilleAbscissae >= Xi(1) - eps
    grevilleAbscissae = Xi(1);
end
if grevilleAbscissae > Xi(end) && grevilleAbscissae <= Xi(end) + eps
    grevilleAbscissae = Xi(end);
end
    
end