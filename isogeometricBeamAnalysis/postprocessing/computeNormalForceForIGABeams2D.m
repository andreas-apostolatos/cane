function N = computeNormalForceForIGABeams2D(knotSpanIndex,p,dHat,parameters,analysis,G,dR)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the normal force at the parametric location where the basis
% functions and the derivatives are computed for the postprocessing of the
% isogeometric Bernoulli and Timoshenko beams in 2D.
%
%      Input :
% knotSpanIndex : The knot span index
%             p : The polynmial degree of the curve
%          dHat : The discrete solution vector
%    parameters : The technical parameters of the beam
%      analysis : Beam analysis type :
%                        'Bernoulli' : isogeometric Bernoulli beam analysis
%                       'Timoshenko' : isogeometric Timoshenko beam analysis
%             G : The base vector of the curve
%            dR : Array containing the first and the second derivatives of the
%                 NURBS curve
%
%        Output :
%             N : The normal force at the given parametric location
%
% Function main body :
%
% 0. Read input
%
% 1. Loop over all the contributions in the knot span
%
% 2. Transform the material law into the Cartesian coordinate system
%
% 3. Compute the normal force
%
%% Function main body

%% 0. Read input

% Initializations
epsilonX = 0;
epsilonY = 0;

%% 1. Loop over all the contributions in the knot span
for b = 0:p
    if strcmp(analysis.type,'Bernoulli')
        epsilonX = epsilonX + dR(b+1,2)*dHat(2*(knotSpanIndex-p+b)-1,1);
        epsilonY = epsilonY + dR(b+1,2)*dHat(2*(knotSpanIndex-p+b),1);
    elseif strcmp(analysis.type,'Timoshenko')
        epsilonX = epsilonX + dR(b+1,2)*dHat(3*(knotSpanIndex-p+b)-2,1);
        epsilonY = epsilonY + dR(b+1,2)*dHat(3*(knotSpanIndex-p+b)-1,1);
    end
end

%% 2. Transform the material law into the Cartesian coordinate system
L = norm(G);
EAxial = parameters.EYoung * parameters.A / L^2;

%% 3. Compute the normal force
N = EAxial*(epsilonX*G(1,1) + epsilonY*G(2,1));

end

