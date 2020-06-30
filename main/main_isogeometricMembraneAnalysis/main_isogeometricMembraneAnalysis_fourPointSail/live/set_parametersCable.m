function parameters = set_parametersCable(EYoung, density, areaCS, prestress)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
% 
% Sets the material and geometrical parameters of the boundary cables
%
%     Input :
%    EYoung : Young's modulus
%   density : Material density
%    areaCS : Cross sectional area
% prestress : Prestress [N/m^2]
%
%     Output :
% parameters : Structure containing the following fields,
%               .E : Young's modulus
%             .rho : Material density
%          .areaCS : Cross sectional area
%       .prestress : Prestress containing the field voigtVector
%
%% Function main body

% Young's modulus
parameters.E = EYoung;

% Poisson ratio
parameters.rho = density;

% Thickness of the membrane
parameters.areaCS = areaCS;

% Prestress of the membrane
parameters.prestress = prestress;

end