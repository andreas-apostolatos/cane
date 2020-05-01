function f = last_windprofil(p, i, xi, Xi, q, j, eta, Eta, CP)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Windprofil for Zone 2, z>4m. qref=0.39N/m^2 ->m^2!!
% Returns Neumann boundary condition that has been recovered from 
% analytical solution
%
%  Input : 
%    p,q : Polynomial degrees
%    i,j : Knot span indeces
% xi,eta : Parametric coordinates on the surface
% Xi,Eta : Knot vectors of the surface
%     CP : The set of control points and weights
%
% Output :
%      f : The value of the analytical load at the surface location (u,v) 
%
%% Function main body

Nu=Bspline_basis_functions(i,p,xi,Xi);
Nv=Bspline_basis_functions(j,q,eta,Eta);
SumNw = 0;
for c = 0:q
    for b = 0:p
        SumNw = SumNw + Nu(b + 1)*Nv(c + 1)*CP(i - p + b, j - q + c, 4);
    end
end

y = 0;
for c = 0:q 
    for b = 0:p
        y = Nu(b + 1)*Nv(c + 1)*CP(i - p + b, j - q + c, 4)*CP(i - p + b, j - q + c, 2)/SumNw + y; 
    end
end

f = 2.1*0.39*(y/10)^0.24;

end
