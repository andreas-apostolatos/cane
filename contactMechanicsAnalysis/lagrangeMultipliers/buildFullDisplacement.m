function [displacement] = buildFullDisplacement(N_unknown,bc,d_reduced)
%% Function documentation
%
% BUILDFULLDISPLACEMENT: All prescribed DoFs are set to be zero. The other
%                        DoFs are assigned to value contained in dred
% 
%              Input :
%          N_unknown : Number of non prescribed DoF
%                 bc : List of indices of prescribed DoFs  
%               dred : Vector with the values for the assignment
%
%             Output :
%       displacement : Vector with the displacement on every DoF
%                     (if this function is called out of an expanded system 
%                     the vector 'displacement' contains also Lagrange multipliers 'dexp')
%
%%
% Dimension of the complete displacement vector
displacement = zeros(N_unknown,1);
% Initialize counters
i=1;
k=1;
for l = 1:N_unknown
    % if we are in a Dirichlet boundary condition location add 0
    if (i<=length(bc) && l==bc(i))
        displacement(l,1) = 0;
        % update counter
        i=i+1;
    % if not add the Control Point displacement 
    else
        displacement(l,1) = d_reduced(k);
        % update counter
        k=k+1;
    end
end

end
 
