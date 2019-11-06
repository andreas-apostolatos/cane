function f = last_windprofil(p,i,u,U,q,j,v,V,CP)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% windprofil f�r Zone 2, z>4m. qref=0.39N/m^2 ->m^2!!
% Returns Neumann boundary condition that has been recovered from 
% analytical solution
%
%  Input : 
%    p,q : Polynomial degrees
%    i,j : Knot span indeces
%    u,v : Parametric coordinates on the surface
%     CP : The set of control points and weights
%
% Output :
%      f : The value of the analytical load at the surface location (u,v) 

%% Function main body

Nu=Bspline_basis_functions(i,p,u,U);
Nv=Bspline_basis_functions(j,q,v,V);
SumNw = 0;
for c = 0:q
  for b = 0:p
    SumNw = Nu(b+1)*Nv(c+1)*CP(i-p+b,j-q+c,4)+SumNw;
  end
end

y = 0;
for c = 0:q 
  for b = 0:p
    y = Nu(b+1)*Nv(c+1)*CP(i-p+b,j-q+c,4)*CP(i-p+b,j-q+c,2)/SumNw+y; 
  end
end

f = 2.1*0.39*(y/10)^0.24;

end
