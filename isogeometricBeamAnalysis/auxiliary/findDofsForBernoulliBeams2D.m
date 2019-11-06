function rb = findDofsForBernoulliBeams2D(rb,xi,direction,CP)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the DOFs which are specified within the chosen parametric
% extension and direction for the isogeometric Benroulli beam element in
% 2D.
%
%     Input :
%        rb : previous set of supports
%        xi : region to be supported (e.g. u=[0 1])
% direction : direction: 1-x, 2-y
%
%    Output :
%        rb : Vector containing the found DOFs
%
% Function layout :
%
% 0. Read input
%
% 1. Iterate and add new DOFs into the array preserving the old ones
%
% 2. Loop over the DOF array and delete double entries
%
%% Function main body

%% 0. Read input

% counter
r = length(rb) + 1;

% number of control points
nxi = length(CP(:,1));

%% 1. Iterate and add new DOFs into the array preserving the old ones
for i = xi(1)*(nxi-1)+1:xi(2)*(nxi-1)+1
    rb(r) = 2*(i-1) + direction;

    % Round to nearest integer
    rb(r) = round(rb(r));

    % Update counter
    r = r + 1;
end

% sort rb and delete double entries
rb = sort(rb);

% Initialize counter
i = 1;

%% 2. Loop over the DOF array and delete double entries
while i < length(rb)
    if rb(i) == rb(i+1)
        rb(i+1)=[];  
        
        % Decrease counter
        i = i - 1;  
    end
    % Increase counter
    i = i + 1;
end

end
