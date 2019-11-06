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
function plot_torus3d(center,Radius,radius,thetaInterval,phiInterval,...
    azimuthallySymmetryAxis,prestress,compPrestress,color)
%% Function documentation
%
% Plots a torus in 3D.
%
%                   Input :
%                  center : The center of the torus in the Cartesian 
%                           coordinates system
%                  Radius : The radius of the torus
%                  radius : The radius of the torus tube
%           thetaInterval : The interval of the angle across the torus
%             phiInterval : The interval of the angle around the torus
% azimuthallySymmetryAxis : Axis of azimuthally symmetry of the torus
%               prestress : Function handle for the computation of the
%                           prestress
%           compPrestress : Component of the prestress to be visualized
%                           (1 or 2)
%                   color : Color for the surface in case no stress is to
%                           be plotted
%
%                  Output :
%                           Graphics
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the theta parametric coordinates
% ->
%    1i. Initialize phi coordinate
%
%   1ii. Loop over all the phi coordinates
%   ->
%        1ii.1. Compute the Cartesian coordinates of the point on the torus
%
%        1ii.2. Update the phi coordinate
%   <-
%
%  1iii. Update the theta coordinate
% <-
%
% 2. Plot the surface
%
%% Function main body

%% 0. Read input

% Evaluation points
noTheta = 50;
noPhi = 50;

% Evaluation increments
dTheta = (thetaInterval(2) - thetaInterval(1))/(noTheta - 1);
dPhi = (phiInterval(2) - phiInterval(1))/(noPhi - 1);

% Function handles for the computation of the Cartesian coordinates
if strcmp(azimuthallySymmetryAxis,'x')
    zeta1 = @(theta,phi) radius*sin(phi);
    zeta2 = @(theta,phi) (Radius + radius*cos(phi))*cos(theta);
    zeta3 = @(theta,phi) (Radius + radius*cos(phi))*sin(theta);
elseif strcmp(azimuthallySymmetryAxis,'y')
    zeta1 = @(theta,phi) (Radius + radius*cos(phi - pi/2))*cos(theta);
    zeta2 = @(theta,phi) radius*sin(phi - pi/2);
    zeta3 = @(theta,phi) (Radius + radius*cos(phi - pi/2))*sin(theta);
elseif strcmp(azimuthallySymmetryAxis,'z')
    zeta1 = @(theta,phi) (Radius + radius*cos(phi))*cos(theta);
    zeta2 = @(theta,phi) (Radius + radius*cos(phi))*sin(theta);
    zeta3 = @(theta,phi) radius*sin(phi);
end

% Initialize array of the Cartesian location of the points on the torus
Torus = zeros(noTheta,noPhi,3);

% Check if prestress value is to be visualized
isPrestress = false;
if ~ischar(compPrestress)
    if ~isempty(compPrestress)
        prestressValues = zeros(noTheta,noPhi,1);
        isPrestress = true;
    end
end

% Initialize theta coordinate
theta = 0;

%% 1. Loop over all the theta parametric coordinates
for iTheta = 1:noTheta
    %% 1i. Initialize phi coordinate
    phi = phiInterval(1);
    
    %% 1ii. Loop over all the phi coordinates
    for iPhi = 1:noPhi
        %% 1ii.1. Compute the Cartesian coordinates of the point on the torus
        Torus(iTheta,iPhi,1) = center(1,1) + zeta1(theta,phi);
        Torus(iTheta,iPhi,2) = center(2,1) + zeta2(theta,phi);
        Torus(iTheta,iPhi,3) = center(3,1) + zeta3(theta,phi);
        
        %% Compute the prestress values
        if isPrestress
            prestressValue = prestress.voigtVector([theta; phi]);
            if compPrestress == 1
                prestressValues(iTheta,iPhi,1) = prestressValue(1,1);
            elseif compPrestress == 2
                prestressValues(iTheta,iPhi,1) = prestressValue(2,1);
            else
                error('The prestress component can either be 1 or 2');
            end
        end
        
        %% 1ii.2. Update the phi coordinate
        phi = phi + dPhi;
    end
    
    %% 1iii. Update the theta coordinate
    theta = theta + dTheta;
end

%% 2. Plot the surface
if isPrestress
    surf(Torus(:,:,1),Torus(:,:,2),Torus(:,:,3),prestressValues,'EdgeColor','none');
    colorbar;
else
    surf(Torus(:,:,1),Torus(:,:,2),Torus(:,:,3),'FaceColor',color,'EdgeColor','none');
end

end