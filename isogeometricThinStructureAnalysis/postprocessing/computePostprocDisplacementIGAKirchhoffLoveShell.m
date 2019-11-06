function d = computePostprocDisplacementIGAKirchhoffLoveShell(p,q,R,dHatActual)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the displacement field d = [dx dy dz]' at the given parametric
% location given the actual Control Point displacement field which affects
% the current knot span for the problem of an isogeometric Kirchhoff-Love
% shell
%
%      Input :
%        p,q : The polynomial degrees in xi- and eta-directions
%          R : The basis functions which affect the current knot span
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

% Initialize counters
counter = 1;
counterDisp = 1;

%% 1. Compute the displacement components iteratively by summing all the contriutions from the basis functions
for c = 0:q
    for b = 0:p
        % Compute the x-component of the displacement field
        d(1,1) = d(1,1) + R(counter,1)*dHatActual(counterDisp,1);

        % Compute the y-component of the displacement field
        d(2,1) = d(2,1) + R(counter,1)*dHatActual(counterDisp + 1,1);
        
        % Compute the z-component of the displacement field
        d(3,1) = d(3,1) + R(counter,1)*dHatActual(counterDisp + 2,1);

        % Update counters
        counter = counter + 1;
        counterDisp = counterDisp + 3; 
    end
end

end

