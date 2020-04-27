function i = findKnotSpan(xi, Xi, numCPs_xi)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% returns the knot span where u lies
% Note: due to rounding errors in u findspan can be ambiguous at knots.
% There's no problem for the inner knots but to get the last knot (special
% case) rounding range must be considered
%
%     Input :
%        xi : The coordinate in the unit interval
%        Xi : The knot vector of the NURBS curve
% numCPs_xi : n = m - p - 1
%
%    Output :
%         i : The knot span where u lies
%
%% Function main body

% Number of knots
m = length(Xi);

% Tolerance
eps = 1e-7;

if norm(xi) < Xi(1) - 1e-7
    xi = Xi(1);
end

% special case: last knot (open knot vector assumed)
if abs(xi-Xi(numCPs_xi+1)) < eps
    i = numCPs_xi;
    return
end

for i = 1:m-1
    if xi < Xi(i+1)
        return
    end
end  

% If no return has been called issue an error
error('xi outside of Xi!');

end
