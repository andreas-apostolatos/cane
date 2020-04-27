function Q = computeShearForceForIGABeams2D ...
    (knotSpanIndex, p, dHat, parameters, analysis, G, dR)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the shear force at the parametric location where the basis
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
%             Q : The shear force at the given parametric location
%
% Function main body :
%
% 0. Read input
%
% 1. Loop over all the contributions in the knot span
%
% 2. Transform the material law into the Cartesian coordinate system
%
% 3. Compute the shear force
%
%% Function main body

%% 0. Read input

% Initializations
if strcmp(analysis.type,'Bernoulli')
    error('The shear force computation for the Bernoulli beam has not yet been implemented');
elseif strcmp(analysis.type,'Timoshenko')
    gammaX = 0;
    gammaY = 0;
    gammaBeta = 0;
end

% Compute the length of the base vector at the parametric location
L = norm(G);

%% 1. Loop over all the contributions in the knot span
for b = 0:p
    if strcmp(analysis.type,'Bernoulli')
        error('The shear force computation for the Bernoulli beam has not yet been implemented');
    elseif strcmp(analysis.type,'Timoshenko')
        gammaX = gammaX + (-G(2,1)/L)*dR(b+1,2)*dHat(3*(knotSpanIndex-p+b)-2,1);
        gammaY = gammaY + (G(1,1)/L)*dR(b+1,2)*dHat(3*(knotSpanIndex-p+b)-1,1);
        gammaBeta = gammaBeta + (-L)*dR(b+1,1)*dHat(3*(knotSpanIndex-p+b),1);
    end
end

%% 2. Transform the material law into the Cartesian coordinate system
EShear = parameters.GShear * parameters.A / L;

%% 3. Compute the shear force
if strcmp(analysis.type,'Bernoulli')
    error('The shear force computation for the Bernoulli beam has not yet been implemented')
elseif strcmp(analysis.type,'Timoshenko')
    Q = EShear*(gammaX + gammaY + gammaBeta);
end

end
