function haveTheSameOrientation = findSubdomainInterfaceOrientation...
    (p1,Xi1,q1,Eta1,CP1,isNURBS1,xicoup1,etacoup1,p2,Xi2,q2,Eta2,CP2,...
    isNURBS2,xicoup2,etacoup2)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Determines whether the coupling surfaces are oriented in the same 
% direction or not in the pshysical space.
%
%                  Input :
%                  p1,p2 : the polynomial degrees in xi-direction for both 
%                          patches
%                  q1,q2 : the polynomial degrees in eta-direction for both 
%                          patches
%                Xi1,Xi2 : the knot vectors in xi-direction for both 
%                          patches
%              Eta1,Eta2 : the knot vectors in eta-direction for both 
%                          patches
%                CP1,CP2 : The control point coordiantes for both patches
%      isNURBS1,isNURBS2 : Flags on whether the bases is a a B-Spline or a 
%                          NURBS 
%        xicoup1,xicoup2 : The coupling region in xi-direction for both 
%                          patches
%      etacoup1,etacoup2 : The coupling region in eta-direction for both 
%                          patches
%     
%                 Output :
% haveTheSameOrientation : Flag determining whether the two coupling curves 
%                          are oriented in the same direction:
%                          1 : they are oriented in the same direction
%                          0 : they are oriented in opposite directions
%
% Function layout :
%
% 0. Read input
%
% 1. Compute three Cartesian locations on the interfaces for bo patches; Start, middle and end points in the parameter space
%
% 2. Decide on whether the interface parametrizations have the same orientation judging by the three points
%
%% Function main body

%% 0. Read input

% Assign a tolerance value
tolerance = 1e1;

% Number of Control Points at xi-,eta- direction for both patches

% 1st patch :
% ___________

nxi1 = length(CP1(:,1,1));
neta1 = length(CP1(1,:,1));

% 2nd patch :
% ___________

nxi2 = length(CP2(:,1,1));
neta2 = length(CP2(1,:,1));
   
%% 1. Compute three Cartesian locations on the interfaces for bo patches; Start, middle and end points in the parameter space

% 1st patch :
% ___________

% Compute the Cartesian coordinates of the start point in the parameter
% space
xi1 = xicoup1(1);
eta1 = etacoup1(1);
xiSpan1 = findKnotSpan(xi1,Xi1,nxi1);
etaSpan1 = findKnotSpan(eta1,Eta1,neta1);
R1 = computeIGABasisFunctionsAndDerivativesForSurface...
    (xiSpan1,p1,xi1,Xi1,etaSpan1,q1,eta1,Eta1,CP1,isNURBS1,0);
xStart1 = computeCartesianCoordinatesOfAPointOnBSplineSurface...
    (xiSpan1,p1,xi1,Xi1,etaSpan1,q1,eta1,Eta1,CP1,R1);

% Compute the Cartesian coordinates of the middle point in the parameter
% space
if etacoup1(1) == etacoup1(2)
    xi1 = (xicoup1(1) + xicoup1(2))/2;
    eta1 = etacoup1(1);
else
    xi1 = xicoup1(1);
    eta1 = (etacoup1(1) + etacoup1(2))/2;
end
xiSpan1 = findKnotSpan(xi1,Xi1,nxi1);
etaSpan1 = findKnotSpan(eta1,Eta1,neta1);
R1 = computeIGABasisFunctionsAndDerivativesForSurface...
    (xiSpan1,p1,xi1,Xi1,etaSpan1,q1,eta1,Eta1,CP1,isNURBS1,0);
xMiddle1 = computeCartesianCoordinatesOfAPointOnBSplineSurface...
    (xiSpan1,p1,xi1,Xi1,etaSpan1,q1,eta1,Eta1,CP1,R1);

% Compute the Cartesian coordinates of the end point in the parameter space
if etacoup1(1) == etacoup1(2)
    xi1 = xicoup1(2);
    eta1 = etacoup1(1);
else
    xi1 = xicoup1(1);
    eta1 = etacoup1(2);
end
xiSpan1 = findKnotSpan(xi1,Xi1,nxi1);
etaSpan1 = findKnotSpan(eta1,Eta1,neta1);
R1 = computeIGABasisFunctionsAndDerivativesForSurface...
    (xiSpan1,p1,xi1,Xi1,etaSpan1,q1,eta1,Eta1,CP1,isNURBS1,0);
xEnd1 = computeCartesianCoordinatesOfAPointOnBSplineSurface...
    (xiSpan1,p1,xi1,Xi1,etaSpan1,q1,eta1,Eta1,CP1,R1);

% 2nd patch :
% ___________

% Compute the Cartesian coordinates of the start point in the parameter
% space
xi2 = xicoup2(1);
eta2 = etacoup2(1);
xiSpan2 = findKnotSpan(xi2,Xi2,nxi2);
etaSpan2 = findKnotSpan(eta2,Eta2,neta2);
R2 = computeIGABasisFunctionsAndDerivativesForSurface...
    (xiSpan2,p2,xi2,Xi2,etaSpan2,q2,eta2,Eta2,CP2,isNURBS2,0);
xStart2 = computeCartesianCoordinatesOfAPointOnBSplineSurface...
    (xiSpan2,p2,xi2,Xi2,etaSpan2,q2,eta2,Eta2,CP2,R2);

% Compute the Cartesian coordinates of the middle point in the parameter
% space
if etacoup2(1) == etacoup2(2)
    xi2 = (xicoup2(1) + xicoup2(2))/2;
    eta2 = etacoup2(1);
else
    xi2 = xicoup2(1);
    eta2 = (etacoup2(1) + etacoup2(2))/2;
end
xiSpan2 = findKnotSpan(xi2,Xi2,nxi2);
etaSpan2 = findKnotSpan(eta2,Eta2,neta2);
R2 = computeIGABasisFunctionsAndDerivativesForSurface...
    (xiSpan2,p2,xi2,Xi2,etaSpan2,q2,eta2,Eta2,CP2,isNURBS2,0);
xMiddle2 = computeCartesianCoordinatesOfAPointOnBSplineSurface...
    (xiSpan2,p2,xi2,Xi2,etaSpan2,q2,eta2,Eta2,CP2,R2);

% Compute the Cartesian coordinates of the end point in the parameter space
if etacoup2(1) == etacoup2(2)
    xi2 = xicoup2(2);
    eta2 = etacoup2(1);
else
    xi2 = xicoup2(1);
    eta2 = etacoup2(2);
end
xiSpan2 = findKnotSpan(xi2,Xi2,nxi2);
etaSpan2 = findKnotSpan(eta2,Eta2,neta2);
R2 = computeIGABasisFunctionsAndDerivativesForSurface...
    (xiSpan2,p2,xi2,Xi2,etaSpan2,q2,eta2,Eta2,CP2,isNURBS2,0);
xEnd2 = computeCartesianCoordinatesOfAPointOnBSplineSurface...
    (xiSpan2,p2,xi2,Xi2,etaSpan2,q2,eta2,Eta2,CP2,R2);

%% 2. Decide on whether the interface parametrizations have the same orientation judging by the three points

if norm(xMiddle1 - xMiddle2) > tolerance
    error('The middle point of the interface parametrizations is not the same');
end
if norm(xStart1 - xStart2) < tolerance && norm(xEnd1 - xEnd2) < tolerance
    haveTheSameOrientation = true;
elseif norm(xStart1 - xEnd2) < tolerance && norm(xEnd1 - xStart2) < tolerance
    haveTheSameOrientation = false;
else
    error('The interface parametrization do not match');
end
