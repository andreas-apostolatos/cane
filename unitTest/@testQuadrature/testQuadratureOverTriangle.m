function testQuadratureOverTriangle(testCase)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
%% Function documentation
%
% Tests the Gauss quadrature over a triangle which is constructed using the
% degenerated quadrilateral for the function,
% f(zeta1,zeta2) = ln(zeta1 + 1)*sin(pi*zeta2/2)
%
% Function layout :
%
% 1. Loop over all quadrature rules
% ->
%    1i. Initialize the linear combination of the Gauss point coordinates with their weights
%
%   1ii. Issue the Gauss point coordinates and weights
%
%  1iii. Loop over all Gauss points and compute the integral
%
%   1iv. Verify the solution
% <-
%
%% Function main body

%% 0. Read input

% Expected solution
expSol = 0.050909793169365;

% Tolerance
tol = 1e-1;
absTol = 1e-15;

%% 1. Loop over all polynomial orders for the Gauss integration over a triangle using the symmetric rule
for iPolOrder = 1:8
    %% 1i. Initialize the linear combination of the Gauss point coordinates with their weights
    integral = 0;
    
    %% 1ii. Issue the Gauss point coordinates and weights
    [GP,GW] = getGaussRuleOnCanonicalTriangle(iPolOrder);
    
    %% 1iii. Loop over all Gauss points and compute the integral
    for iGP = 1:length(GP(:,1))
        integrand = log(GP(iGP,1) + 1)*sin(pi*GP(iGP,2)/2);
        integral = integral + integrand*GW(iGP,1);
    end
    
    %% 1iv. Verify the solution
    tolStricktened = tol*10^(-iPolOrder + 2);
    testCase.verifyEqual(integral,expSol,'AbsTol',tolStricktened);
end

%% 2. Loop over all the quadrature rules for the Gauss integration over a triangle using the degenerated quadrilateral
for iQuadrature = 1:10
    %% 2i. Initialize the linear combination of the Gauss point coordinates with their weights
    integral = 0;
    
    %% 2ii. Issue the Gauss point coordinates and weights
    noGP = iQuadrature^2;
    [GP,GW] = getGaussRuleOnCanonicalTriangleWithDegeneratedQuadrilateral(noGP);
    
    %% 2iii. Loop over all Gauss points and compute the integral
    for iGP = 1:length(GP(:,1))
        integrand = log(GP(iGP,1) + 1)*sin(pi*GP(iGP,2)/2);
        integral = integral + integrand*GW(iGP,1);
    end
    
    %% 2iv. Verify the solution
    tolStricktened = tol*10^(-3*iQuadrature/2 + (iQuadrature == 1)*(2/2));
    testCase.verifyEqual(integral,expSol,'AbsTol',tolStricktened*(tolStricktened > absTol) + absTol*(tolStricktened <= absTol));
end

end
