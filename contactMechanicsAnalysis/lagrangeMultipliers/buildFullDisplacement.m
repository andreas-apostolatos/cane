function displacement = buildFullDisplacement(nDOF,homDBC,displacement_red)
%% Function documentation
%
% BUILDFULLDISPLACEMENT: All prescribed DoFs are set to be zero. The other
%                        DoFs are assigned to value contained in dred
% 
%              Input :
%               nDOF : Number of non prescribed DoF
%             homDBC : List of indices of prescribed DoFs  
%   displacement_red : Vector with the values for the assignment
%
%             Output :
%       displacement : Vector with the displacement on every DoF
%                     (if this function is called out of an expanded system 
%                     the vector displacement contains also Lagrange multipliers 'dexp')
%
%%
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