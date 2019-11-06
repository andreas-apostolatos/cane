function [d,gradd,gradgradd] = ...
    computePostprocDisplacementAndGradientIGAKirchhoffLoveShell...
    (p,q,dR,dHatActual)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the displacement field d = [dx dy dz]', and its first and second
% order gradient namely, 
%
%       grad(d) = | d1,1 d1,2 |
%                 | d2,1 d2,2 |
%
% grad(grad(d)) = | d1,11 d1,21 d1,12 d1,22 |
%                 | d2,11 d2,21 d2,12 d2,22 |
% 
% given the discrete displacement vector that affects the current knot
% span.
%
%      Input :
%        p,q : The polynomial degrees in xi- and eta-directions
%         dR : The basis functions which affect the current knot span and
%              their first and second derivatives, namely,
%              dR = [R dR/dxi d^2R/dxi dR/deta d^2R/deta/dxi d^2R/deta]
% dHatActual : The Control Point displacement field affecting the current
%              knot span
%
%     Outout :
%          d : the displacement field d = [dx dy dz]' at the given parametric
%              location
%
% Function layout :
%
% 0. Read input
%
% 1. Compute the displacement components iteratively by summing all the contriutions from the basis functions
%
%% Function main body

%% 0. Read input

% Initialize output array
d = zeros(3,1);
gradd = zeros(2,2);
gradgradd = zeros(2,4);

% Initialize counters
counter = 1;
counterDisp = 1;

%% 1. Compute the displacement components iteratively by summing all the contriutions from the basis functions
for c = 0:q
    for b = 0:p
        %% 1i. Compute the displacement field itself
        % Compute the x-component of the displacement field
        d(1,1) = d(1,1) + R(counter,1)*dHatActual(counterDisp,1);

        % Compute the y-component of the displacement field
        d(2,1) = d(2,1) + R(counter,1)*dHatActual(counterDisp + 1,1);
        
        % Compute the z-component of the displacement field
        d(3,1) = d(3,1) + R(counter,1)*dHatActual(counterDisp + 2,1);
        
        %% 1ii. Compute the gradient of the displacement field
        
        %% 1iii. Compute the second order gradient of the displacement field

        % Update counters
        counter = counter + 1;
        counterDisp = counterDisp + 3; 
    end
end

end

