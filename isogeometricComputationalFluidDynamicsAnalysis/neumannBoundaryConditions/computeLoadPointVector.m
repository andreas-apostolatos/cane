function fl = computeLoadPointVector ...
    (fl_old, xi, eta, p, q, Xi, Eta, CP, f, dir)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% returns the consistent nodal forces to a point load f.
%
%     Input :
%    fl_old : Existing force vector
%    xi,eta : Parametric extension of the location where the load is 
%             applied
%    Xi,Eta : Knot vectors along the xi- and eta-directions
%        CP : Control point coordinates and weights
%         f : Load magnitude
%       dir : Load direction,
%               1 : x, 
%               2 : y
%
%    Output :
%        fl : Updated force vector
%
%% Function main body

numCPs_xi = length(CP(:, 1, 1));
numCPs_eta = length(CP(1, :, 1));

i = findspan(xi, Xi, numCPs_xi);
j = findspan(eta, Eta, numCPs_eta);

% Compute the NURBS basis functions in 2D
R = nurbs_basis_functions2D(i, p, xi, Xi, j, q, eta, Eta, CP);

% Initialize the load vector
Fl = zeros(numCPs_xi, numCPs_eta, 2);

% Compute the load as a linear combination of the basis functions at the
% inspection point

% initialize counter
k = 1;

for c = 0:q
    for b = 0:p
        
        Fl(i-p+b,j-q+c,dir) = R(k)*f;
        
        % Update counter
        k = k + 1;
    end
end

fl = make_fl_dof(Fl);
if isvector(fl_old);  fl = fl + fl_old;   end

end
