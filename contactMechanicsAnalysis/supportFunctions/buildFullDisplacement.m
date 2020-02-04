function displacement = buildFullDisplacement(nDOF,homDBC,displacement_red)
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
%              nDOF : Number of non-prescribed DOFs
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
displacement = zeros(nDOF,1);

% Initialize counters
i=1;
k=1;

% loop over degrees of freedom
for l = 1:nDOF
    % if we are in a Dirichlet boundary condition location add 0
    if (i<=length(homDBC) && l==homDBC(i))
        displacement(l,1) = 0;
        % update counter
        i=i+1;
    % if not add the Control Point displacement 
    else
        displacement(l,1) = displacement_red(k);
        % update counter
        k=k+1;
    end
end

end