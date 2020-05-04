function testBSplineUtilityFunctions(testCase)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Tests the B-Spline utility functions
%
% - Projection on a quarter of a sphere geometry on a regular and on an
%   irregular parametric location
%
% Function layout :
%
% 0. Read input
%
% 1. Define the CAD geometry in terms of NURBS
%
% 2. p-refine the surface
%
% 3. h-refine the surface
%
% 4. Project the points on the surfaces
%
% 5. Verify the results
%
%% Function main body

%% 0. Read input

% Absolute tolerances
absTol = 1e-15;
absTolRel9 = absTol*1e9;
absTolRel12 = absTol*1e12;

%% 1. Define the CAD geometry in terms of NURBS

% Global variables
Radius = 0.075;
        
% Polynomial degrees
p1 = 2;
q1 = 2;

% Knot vectors
Xi1 = [0 0 0 1 1 1];
Eta1 = [0 0 0 1 1 1];

% Control Point coordinates

% x-coordinates
CP1(:,:,1) = [0      0      0
              Radius Radius Radius
              Radius Radius Radius];

% y-coordinates
CP1(:,:,2) = [-Radius -Radius 0
              -Radius -Radius 0
              0       0       0];

% z-coordinates for an exact semisphere
CP1(:,:,3) = [0 Radius Radius
              0 Radius Radius
              0 0      0];

% Weights
weight = sqrt(2)/2;
CP1(:,:,4) = [1      weight  1
              weight weight^2 weight 
              1      weight  1];

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBS1 = 0;
nxi1 = length(CP1(:,1,1));
neta1 = length(CP1(1,:,1));
for i = 1:nxi1
    for j = 1:neta1
        if CP1(i,j,4)~=1
            isNURBS1 = 1;
            break;
        end
    end
    if isNURBS1
        break;
    end
end

%% 2. p-refine the surface

tp = 1;
tq = 2;
[Xi1,Eta1,CP1,p1,q1] = degreeElevateBSplineSurface(p1,q1,Xi1,Eta1,CP1,tp,tq,'');

%% 3. h-refine the surface

noXi = 3;
noEta = 2;
[Xi1,Eta1,CP1] = knotRefineUniformlyBSplineSurface(p1,Xi1,q1,Eta1,CP1,noXi,noEta,'');

%% 4. Project the points on the surfaces

% Set an initial guess for both projections
xi0 = .5;
eta0 = .5;

% Point to be projected at a regular parametric location
P_reg = [0.06
         -0.06
         0.06];
     
% Set the projection parameters
propNewtonRaphson.maxIt = 40;
propNewtonRaphson.eps = 1e-12;
     
% Closest point projection
[xi_reg,eta_reg,Projected_reg,isProjected_reg,noIter_reg] = ...
    computeNearestPointProjectionOnBSplineSurface...
    (P_reg,p1,Xi1,q1,Eta1,CP1,isNURBS1,xi0,eta0,propNewtonRaphson);

% Test the projection at a regular parametric location
expSol_xi_reg = 0.397197676620219;
expSol_eta_reg = 0.500000000000000;
expSol_Projected_reg = [0.043301270189222
                        -0.043301270189222
                        0.043301270189222];

% Point to be projected at an irregular parametric location
P_irreg = [0.1; -0.00005; 0.00005];

% Set the projection parameters
propNewtonRaphson.maxIt = 100;
propNewtonRaphson.eps = 1e-12;

% Closest point projection
[xi_irreg,eta_irreg,Projected_irreg,isProjected_irreg,noIter_irreg] = ...
    computeNearestPointProjectionOnBSplineSurface...
    (P_irreg,p1,Xi1,q1,Eta1,CP1,isNURBS1,xi0,eta0,propNewtonRaphson);

% Test the projection at a regular parametric location
expSol_xi_irreg = 1.000000000000000;
expSol_eta_irreg = 0.500000000000001;
expSol_Projected_irreg = [0.075000000000000
                          -0.000000000000000
                          0.000000000000000];

%% 5. Verify the results
                      
expSol_isProjected_reg = true;
expSol_noIter_reg = 4;
testCase.verifyEqual(xi_reg,expSol_xi_reg,'AbsTol',absTol);
testCase.verifyEqual(eta_reg,expSol_eta_reg,'AbsTol',absTol);
testCase.verifyEqual(Projected_reg,expSol_Projected_reg,'AbsTol',absTol);
testCase.verifyEqual(isProjected_reg,expSol_isProjected_reg,'AbsTol',absTol);
testCase.verifyEqual(isProjected_reg,expSol_isProjected_reg,'AbsTol',absTol);
testCase.verifyEqual(noIter_reg,expSol_noIter_reg,'AbsTol',absTol);
                      
expSol_isProjected_irreg = true;
expSol_noIter_irreg = 4;
testCase.verifyEqual(xi_irreg,expSol_xi_irreg,'AbsTol',absTolRel12);
testCase.verifyEqual(eta_irreg,expSol_eta_irreg,'AbsTol',absTolRel9);
testCase.verifyEqual(Projected_irreg,expSol_Projected_irreg,'AbsTol',absTolRel12);
testCase.verifyEqual(isProjected_irreg,expSol_isProjected_irreg,'AbsTol',absTol);
testCase.verifyEqual(isProjected_irreg,expSol_isProjected_irreg,'AbsTol',absTol);
testCase.verifyEqual(noIter_irreg,expSol_noIter_irreg,'AbsTol',absTol);

end
