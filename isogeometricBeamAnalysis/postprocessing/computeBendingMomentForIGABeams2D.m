function M = computeBendingMomentForIGABeams2D ...
    (knotSpanIndex, p, dHat, parameters, analysis, G, dG, dR)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the bending moment at the parametric location where the basis
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
%            dG : The acceleration vector
%            dR : Array containing the first and the second derivatives of the
%                 NURBS curve
%
%        Output :
%             M : The bending moment at the given parametric location
%
% Function main body :
%
% 0. Read input
%
% 1. Loop over all the contributions in the knot span
%
% 2. Transform the material law into the Cartesian coordinate system
%
% 3. Compute the bending moment
%
%% Function main body

%% 0. Read input

% Initializations
if strcmp(analysis.type,'Bernoulli')
    kappaX = 0;
    kappaY = 0;
elseif strcmp(analysis.type,'Timoshenko')
	kappaX = 0;
    kappaY = 0;
    kappaBeta = 0;
end

% Compute the length of the base vector at the parametric location
L = norm(G);

%% 1. Loop over all the contributions in the knot span
if strcmp(analysis.type,'Bernoulli')
    % Necessary variables for the curvature
    d1 = (1/L^2)*G(1,1)*G(2,1);
    d2 = (1/L^2)*G(2,1)^2-1;
    d3 = 1-(1/L^2)*G(1,1)^2;

    % Variables alpha's
    alpha1 = (dG(1,1)*d1+dG(2,1)*d3)/L ;
    alpha2 = (dG(1,1)*d2-dG(2,1)*d1)/L ;

    % Loop over all CP's affecting the current knot span
    for b = 0:p
        if strcmp(analysis.type,'Bernoulli')
            kappaX = kappaX + (-alpha1*dR(b+1,2)+G(2,1)*dR(b+1,3)/L)*dHat(2*(knotSpanIndex-p+b)-1,1);
            kappaY = kappaY + (-alpha2*dR(b+1,2)-G(1,1)*dR(b+1,3)/L)*dHat(2*(knotSpanIndex-p+b),1);
        elseif strcmp(analysis.type,'Timoshenko')
            kappaX = kappaX + (-alpha1*dR(b+1,2)+G(2,1)*dR(b+1,3)/L)*dHat(3*(knotSpanIndex-p+b)-2,1);
            kappaY = kappaY + (-alpha2*dR(b+1,2)-G(1,1)*dR(b+1,3)/L)*dHat(3*(knotSpanIndex-p+b)-1,1);
        end
    end
elseif strcmp(analysis.type,'Timoshenko')
    % Variable alpha
    alpha = (1/L^3)*(dG(1)*G(1,1)+dG(2)*G(2,1));

    % Differentiation of A2 with respect to theta^1
    A2xKomma1 = - (1/L)*dG(2)+alpha*G(2,1); 
    A2yKomma1 = (1/L)*dG(1)-alpha*G(1,1);

    % Variable beta
    beta = - G(1,1) * A2yKomma1 + G(2,1) * A2xKomma1 ;

    % Loop over all CP's affecting the current knot span
    for b = 0:p
        kappaX = kappaX + A2xKomma1*dR(b+1,2)*dHat(3*(knotSpanIndex-p+b)-2,1);
        kappaY = kappaY + A2yKomma1*dR(b+1,2)*dHat(3*(knotSpanIndex-p+b)-1,1);
        kappaBeta = kappaBeta + (-L*dR(b+1,2)+beta*dR(b+1,1))*dHat(3*(knotSpanIndex-p+b),1);
    end
end

%% 2. Transform the material law into the Cartesian coordinate system
EBending = parameters.EYoung*parameters.I / L^2;

%% 3. Compute the bending moment
if strcmp(analysis.type,'Bernoulli')
    M = EBending*(kappaX + kappaY);
elseif strcmp(analysis.type,'Timoshenko')
    M = EBending*(kappaX + kappaY + kappaBeta);
end

