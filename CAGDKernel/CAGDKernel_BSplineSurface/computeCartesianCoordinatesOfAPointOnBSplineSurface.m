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
function S = computeCartesianCoordinatesOfAPointOnBSplineSurface...
    (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,R)
%% Function documentation
%
% Returns the X,Y,Z coordinates of the point on the NURBS surface 
% corresponding to xi,eta parametric coordinates in the NURBS parent domain
%
%          Input :
%            p,q : The polynomial degrees in xi-, eta-directions
% xiSpan,etaSpan : The knot span indices in xi-, eta-directions
%         xi,eta : The surface parameters in xi-, eta-directions
%         Xi,Eta : knot vectors in xi-, eta-directions
%             CP : Control Point coordinates and weights
%              R : The NURBS basis functions at (xi,eta)
%
%         Output :
%              S : the Cartesian coordinates of the point with curvilinear 
%                  coordinates (xi,eta)
%
% Function layout
%
% 0. Read input
%
% 1. Loop and summation over all the NURBS basis functions
%
%% Function main body

%% 0. Read input

% Compute the xiSpan, etaSpan knot span indices
if xiSpan == 0
    xiSpan = findKnotSpan(xi,Xi,length(CP(:,1,1)));  
end
if etaSpan == 0
    etaSpan = findKnotSpan(eta,Eta,length(CP(1,:,1)));  
end

% Initialize output array
S = zeros(3,1);

% Initialize iterator
k = 0;

%% 1. Loop and summation over all the NURBS basis functions
for c = 0:q 
    for b = 0:p
        % Update counter
        k = k+1;
        
        % Compute the coordinates iteratively
        S(1,1) = R(k,1)*CP(xiSpan-p+b,etaSpan-q+c,1) + S(1,1);
        S(2,1) = R(k,1)*CP(xiSpan-p+b,etaSpan-q+c,2) + S(2,1);
        S(3,1) = R(k,1)*CP(xiSpan-p+b,etaSpan-q+c,3) + S(3,1); 
    end
end

end