function [n,t] = normal_and_tangent_vectors_to_the_surface_boundary(i,p,u,U,j,q,v,V,CP,is_on_u)
%% Function documentation
%
% Returns the normal and the tangent vectors to the shell boundary, given
% the parametric locations of the boundary location. It should be
% explicitly specified by the last argument of the function in which 
% parametric line the boundary lies on 
%
%   Input : 
%     i,j : Knot span indices
%     p,q : The polynomial degrees in the u,v-directions
%     u,v : The parametric locations on the suface boundary only. Interior
%           locations are not allowed and result into an error in the call 
%           of the function
%     U,V : The knot vectors in u,v-directions
%      CP : The set of Control Points and coordinates
% is_on_u : The 
%
%  Output :
%       n : The normal to the surface boundary vector
%       t : The tangent to the surface vector (always oriented together 
%           with the tanget to the boundary line base vector)
%
% Function layout :
%
% 1. Compute the base vectors at the parametric location
%
% 2. Compute the main normal to surface vector
%
% 3. Check on which boundary the parametric location is located
%
%     3i. Compute the surface normal vector at the parametric location
%
%     3ii. Compute the normalized tangent and normal vector
%
%
%% Function main body

%% 1. Compute the base vectors at the parametric location
[g1,g2] = base_vectors2D(p,i,u,U,q,j,v,V,CP);

%% 2. Compute the main normal to surface vector
A3 = cross_product(g1,g2);
A3 = A3/norm(A3);

%% 3. Check on which boundary the parametric location is located
if is_on_u
    %% ii. Compute the normalized tangent and normal vector
    if v==V(1)
        % Tangent vector
        t = g1/norm(g1);
        % Normal vector
        n = cross_product(t,A3);
    elseif v==V(length(V))
        % Tangent vector
        t = -g1/norm(g1);
        % Normal vector
        n = cross_product(t,A3);
    end
else
    %% ii. Compute the normalized tangent vector
    if u==U(1)
        % Tangent vector
        t = g2/norm(g2);
        % Normal vector
        n = cross_product(A3,t);
    elseif u==U(length(U))
        % Tangent vector
        t = -g2/norm(g2);
        % Normal vector
        n = cross_product(A3,t);
    end
end

end

