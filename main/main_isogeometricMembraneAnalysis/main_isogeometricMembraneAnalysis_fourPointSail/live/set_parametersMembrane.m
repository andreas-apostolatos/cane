function parameters = set_parametersMembrane(EYoung, nue, density, thickness, prestress)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
% 
% Sets the material and geometrical parameters of the membrane
%
%     Input :
%    EYoung : Young's modulus
%       nue : Poisson ratio
%   density : Material density
% thickness : Membrane lateral thickness
% prestress : Prestress containing the following field,
%               .voigtVector: [sigma11; sigma22; sigma12] [N/m^2]
%
%     Output :
% parameters : Structure containing the following fields,
%               .E : Young's modulus
%             .nue : Poisson ratio
%               .t : Lateral thickness
%             .rho : Material density
%       .prestress : Prestress containing the field voigtVector
%
%% Function main body

% Young's modulus
parameters.E = EYoung;

% Poisson ratio
parameters.nue = nue;

% Thickness of the membrane
parameters.t = thickness;

% Density of the membrane
parameters.rho = density;

% Prestress of the membrane
parameters.prestress = prestress;

end