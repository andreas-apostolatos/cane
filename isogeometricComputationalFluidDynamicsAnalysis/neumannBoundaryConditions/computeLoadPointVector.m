%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _______________________________________________________               %
%   _______________________________________________________               %
%                                                                         %
%   Technische Universität München                                        %
%   Lehrstuhl für Statik, Prof. Dr.-Ing. Kai-Uwe Bletzinger               %
%   _______________________________________________________               %
%   _______________________________________________________               %
%                                                                         %
%                                                                         %
%   Authors                                                               %
%   _______________________________________________________________       %
%                                                                         %
%   Dipl.-Math. Andreas Apostolatos    (andreas.apostolatos@tum.de)       %
%   Dr.-Ing. Roland Wüchner            (wuechner@tum.de)                  %
%   Prof. Dr.-Ing. Kai-Uwe Bletzinger  (kub@tum.de)                       %
%   _______________________________________________________________       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fl = computeLoadPointVector(fl_old,u,v,p,q,U,V,CP,f,dir)
%% Function documentation
%
% returns the consistent nodal forces to a point load f.
% Parameters:
%     fl_old      existing force vector
%     fl:         updated force vector
%     u,v:        load position
%     f:          point load
%     dir:        direction of f  1=x, 2=y
%
%% Function main body

nu = length(CP(:,1,1));
nv = length(CP(1,:,1));

i = findspan(u,U,nu);
j = findspan(v,V,nv);

% Compute the NURBS basis functions in 2D
R = nurbs_basis_functions2D(i,p,u,U,j,q,v,V,CP);

% Initialize the load vector
FL = zeros(nu,nv,2);

% Compute the load as a linear combination of the basis functions at the
% inspection point

% initialize counter
k = 1;

for c = 0:q
    for b = 0:p
        
        FL(i-p+b,j-q+c,dir) = R(k)*f;
        
        % Update counter
        k = k + 1;
    end
end

fl = make_fl_dof(FL);
if isvector(fl_old);  fl = fl + fl_old;   end

end