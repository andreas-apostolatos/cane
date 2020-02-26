function displacement = buildFullDisplacement(nDOFsFull,homDBC,displacement_red)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Marko Leskovar
%                  Andreas Apostolatos
%
% Date : 04.02.2020
%
%% Function documentation
%
% All prescribed DOFs are set to be zero. The other DOFs are assigned to
% value contained in displacement_red(uced)
% 
%             Input :
%         nDOFsFull : Number of non-prescribed DOFs
%            homDBC : List of indices of prescribed DOFs  
%  displacement_red : Vector with the values for the assignment
%
%            Output :
%      displacement : Vector with the displacement on every DOF (if this
%                     function is called out of an expanded system the
%                     vector displacement contains also Lagrange multipliers
%
%% Function main body

% Dimension of the complete displacement vector
displacement = zeros(nDOFsFull,1);

% Initialize counters
n=1;
k=1;

% loop over degrees of freedom
for i = 1:nDOFsFull
    % if we are in a Dirichlet boundary condition location add 0
    if (n<=length(homDBC) && i==homDBC(n))
        displacement(i,1) = 0;
        % update counter
        n=n+1;
    % if not add the Control Point displacement 
    else
        displacement(i,1) = displacement_red(k);
        % update counter
        k=k+1;
    end
end

end