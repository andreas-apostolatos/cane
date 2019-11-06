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
function [n,t] = computeNormalAndTangentVectorsToBSplineBoundary...
    (xi,Xi,eta,Eta,GXi,GEta,A3,isOnXi)
%% Function documentation
%
% Returns the normal and the tangent vectors to the shell boundary, given
% the parametric locations of the boundary location given the base vectors. 
% It should be explicitly specified by the last argument of the function in 
% which parametric line the boundary lies on 
%
%    Input : 
%   xi,eta : The parametric locations on the suface boundary only. Interior
%            locations are not allowed and result into an error in the call 
%            of the function
%   Xi,Eta : The knot vectors in u,v-directions
% GXi,GEta : The covariant base vectors on the B-Spline surface
%       A3 : The surface normal base vector
%   isOnXi : The parametric line where the boundary lies on 
%
%   Output :
%        n : The normal to the surface boundary vector
%        t : The tangent to the surface vector (always oriented together 
%            with the tanget to the boundary line base vector)
%
% Function layout :
%
% 1. Compute the tangent to the boundary vector to have always counterclockwise direction
%
% 2. Compute the normal to the boundary vector to point always outwards
%
%% Function main body

%% 1. Compute the tangent to the boundary vector to have always counterclockwise direction
if isOnXi
    if eta == Eta(1)
        t = GXi/norm(GXi);
    elseif eta == Eta(length(Eta))
        t = - GXi/norm(GXi);
    else
        error('Check the coupling or loading extension');
    end
else
    if xi == Xi(1)
        t = - GEta/norm(GEta);
    elseif xi == Xi(length(Xi))
        t = GEta/norm(GEta);
    else
        error('Check the coupling or loading extension');
    end
end

%% 2. Compute the normal to the boundary vector to point always outwards
n = cross(t,A3);

end

