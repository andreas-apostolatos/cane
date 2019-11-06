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
function [propMembrane,propCables] = computeHomogenizedMaterialPropertiesOlympiadach...
    (propPlexiglas,propSteelCables,propMembrane,characteristicLength,noCables)
%% Function documentation
%
% Returns the homogenized material properties for the membrane and the
% boundary cables corresponding to Olympiadach of the stadium in Munich.
%
%                   Input :
%           propPlexiglas : Properties for the plexiglas plates
%                               .EYoung : Young's modulus
%                                  .nue : Poisson's ratio
%                                  .rho : Density
%                            .thickness : Thickness
%         propSteelCables : Properties for the steel cables
%                               .EYoung : Young's modulus
%                                  .nue : Poisson's ratio
%                                  .rho : Density
%                             .diameter : Diameter of the cables' cross section
%                                 .area : Area of the cables' cross section
%            propMembrane : Properties for the membrane model
%                                  .nue : Poisson's ratio
%    characteristicLength : Characteristic length of the element where the
%                           homogenization takes place
%                noCables : Number of cables per characteristic element
%
%                  Output :
%            propMembrane : Updated properties for the membrane model with
%                               .EYoung : Homogenized Young's modulus
%                                  .nue : Poisson's ratio
%                            .thickness : Thickness of the membrane
%              propCables : Update properties for the cables' model with
%                               .EYoung : Updated Youngs' modulus
%                             .diameter : Diameter of the cables' cross section
%                                 .area : Area of the cables' cross section
%
% Function layout :
%
% 1. Compute the homogenized thickness of the membrane
%
% 2. Compute the area and the diameter of the boundary cables
%
% 3. Compute the volume fractions
%
% 4. Compute the homogenized material parameters using the volume fraction method
%
% 5. Compute the Young's modulus of the boundary cables
%
% 6. Compute the density for the membrane
%
% 7. Compute the density for the cables
%
%% Basic material properties

%% 1. Compute the homogenized thickness of the membrane
propMembrane.thickness = propPlexiglas.thickness + 2*propSteelCables.area/characteristicLength;

%% 2. Compute the area and the diameter of the boundary cables
propCables.area = 10*propMembrane.thickness;
propCables.diameter = 2*sqrt(propCables.area/pi);
diameterCablesActual = 5*propSteelCables.diameter;
areaCablesActual = pi*(diameterCablesActual/2)^2;

%% 3. Compute the volume fractions

% Compute the cross sectional area of the plexiglas in the reference square
crossSectionalArea_Plexiglas = characteristicLength*propPlexiglas.thickness;

% Compute the cross sectional area of the cables in the reference square
crossSectionalArea_Cables = (noCables/2)*propSteelCables.area;

% Compute the cross sectional area
crossSectionalArea = crossSectionalArea_Plexiglas + crossSectionalArea_Cables;

% Compute the corresponding volume fractions
volumeFraction_Plexiglas = crossSectionalArea_Plexiglas/crossSectionalArea;
volumeFraction_Cables = crossSectionalArea_Cables/crossSectionalArea;

%% 4. Compute the homogenized material parameters using the volume fraction method

% Homogenized Young's modulus
E_homogenized = propPlexiglas.EYoung*volumeFraction_Plexiglas + propSteelCables.EYoung*volumeFraction_Cables;

% Young modulus for the model
propMembrane.EYoung = E_homogenized/(propMembrane.thickness*characteristicLength)*crossSectionalArea;

%% 5. Compute the Young's modulus of the boundary cables
propCables.EYoung = propSteelCables.EYoung*areaCablesActual/propCables.area;

%% 6. Compute the density for the membrane
massPerSquare75 = propPlexiglas.rho*propPlexiglas.thickness*characteristicLength^2 + ...
    noCables*propSteelCables.rho*propSteelCables.area*characteristicLength;
massPerSquare75 = massPerSquare75 + 0.15*massPerSquare75;
propMembrane.rho = massPerSquare75/characteristicLength^2/propMembrane.thickness;

%% 7. Compute the density for the cables
propCables.rho = propSteelCables.rho*areaCablesActual/propCables.area;

end