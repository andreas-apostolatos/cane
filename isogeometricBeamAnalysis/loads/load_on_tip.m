function F = load_on_tip(CP,amplitude,dof,analysis)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
% Gets the force amplitude and direction on tip and 
% returns the load vector  

% Input :
% CP         : The Control Points of the NURBS curve
% amplitude  : Pressure
% dof        : Horizontal or vertical load
% method     : Bernoulli or Timoshenko technical beam theory 
%
% Output :
% f : Load vector 

% Dof=1 : Horizontal force
% Dof=2 : Vertical force
if  strcmp(analysis.type,'Bernoulli')
    F = zeros(size(CP,1)*2,1) ;
    F(size(CP,1)*2+dof-2,1) = amplitude ; 
elseif strcmp(analysis.type,'Timoshenko')
    F = zeros(size(CP,1)*3,1) ;
    F(size(CP,1)*3+dof-3,1) = amplitude ; 
end


