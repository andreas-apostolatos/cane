function plot_ControlPolygonBSplineCurveOnBSplineSurface ...
    (CP, p_surface, q_surface, Xi_surface, Eta_surface, CP_surface, ...
    isNURBS_surface)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documenation
% 
% Plots Control Points and Control Polygon corresponding to a B-Spline 
% curve
%
%           Input :
%              CP : The control points for a B-Spline curve
%    p*,q_surface : The polynomial degrees of the B-Spline surface
% Xi*,Eta_surface : The knot vectors of the B-Spline surface
%      CP_surface : The Control Point coordinates and the weights of the
%                   B-Spline surface which are living in the physical space
% isNURBS_surface : Flag on whether the surface basis is a B-Spline or a 
%                   NURBS
%
%          Output :
%                   graphics
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the Control Points on the curve
% ->
%    1i. Find the knot spans
%
%   1ii. Compute the NURBS basis functions of the surface
%
%  1iii. Compute the Cartesian coordinates of the curve's Control Points
%
%   1iv. Plot the Control Points
%
%    1v. Initialize the line parameter, the line parameter step size and the segment's left corner
%
%   1vi. Loop over all the segments of the line
%   ->
%        1vi.1. Compute the segment's right corner
%
%        1vi.2. Compute the knot span indices of the left and the right segment's corners
%
%        1vi.3. Compute the surface's basis functions at the left and the right segment's corners
%
%        1vi.4. Compute the Cartesian coordinates of the left and the right segment's corners
%
%        1vi.5. Plot the line segment
%
%        1vi.6. Assign the right segment's corner to the left one
%
%        1vi.7. Update the line parameter
%   <-
% <-
%
%% Function main body

%% 0. Read input

% Number of grid points
noGridPoints = 10;

%% 1. Loop over all the Control Points on the curve
for iCP = 1:length(CP(:, 1)) - 1
    %% 1i. Find the knot spans
    xiSpanPrevious = findKnotSpan(CP(iCP,1), Xi_surface, length(CP_surface(:,1,1)));
    etaSpanPrevious = findKnotSpan(CP(iCP,2), Eta_surface, length(CP_surface(1,:,1)));
    xiSpanNext = findKnotSpan(CP(iCP + 1,1), Xi_surface, length(CP_surface(:,1,1)));
    etaSpanNext = findKnotSpan(CP(iCP + 1,2), Eta_surface, length(CP_surface(1,:,1)));
    
    %% 1ii. Compute the NURBS basis functions of the surface
    RPrevious = computeIGABasisFunctionsAndDerivativesForSurface...
        (xiSpanPrevious,p_surface,CP(iCP,1),Xi_surface,etaSpanPrevious,q_surface,CP(iCP,2),...
        Eta_surface,CP_surface,isNURBS_surface,0);
    RNext = computeIGABasisFunctionsAndDerivativesForSurface...
        (xiSpanNext,p_surface,CP(iCP + 1,1),Xi_surface,etaSpanNext,q_surface,CP(iCP + 1,2),...
        Eta_surface,CP_surface,isNURBS_surface,0);
    
    %% 1iii. Compute the Cartesian coordinates of the curve's Control Points
    SPrevious = computeCartesianCoordinatesOfAPointOnBSplineSurface...
        (xiSpanPrevious,p_surface,CP(iCP,1),Xi_surface,etaSpanPrevious,q_surface,CP(iCP,2),Eta_surface,CP_surface,RPrevious);
    SNext = computeCartesianCoordinatesOfAPointOnBSplineSurface...
        (xiSpanNext,p_surface,CP(iCP + 1,1),Xi_surface,etaSpanNext,q_surface,CP(iCP + 1,2),Eta_surface,CP_surface,RNext);
    
    %% 1iv. Plot the Control Points
    plot3(SPrevious(1,1), SPrevious(2,1), SPrevious(3,1),'--omagenta');
    plot3(SNext(1,1), SNext(2,1), SNext(3,1),'--omagenta');
    
    %% 1v. Initialize the line parameter, the line parameter step size and the segment's left corner
    lambda = 0;
    dlambda = 1/(noGridPoints - 1);
    xiEtaPrevious = CP(iCP,:)';
    
    %% 1vi. Loop over all the segments of the line
    for iGrid = 1:noGridPoints
        %% 1vi.1. Compute the segment's right corner
        xiEtaNext = (1 - lambda)*CP(iCP,1:3)' + lambda*CP(iCP + 1,1:3)';
        
        %% 1vi.2. Compute the knot span indices of the left and the right segment's corners
        xiSpanPrevious = findKnotSpan(xiEtaPrevious(1,1), Xi_surface, length(CP_surface(:,1,1)));
        etaSpanPrevious = findKnotSpan(xiEtaPrevious(2,1), Eta_surface, length(CP_surface(1,:,1)));
        xiSpanNext = findKnotSpan(xiEtaNext(1,1), Xi_surface, length(CP_surface(:,1,1)));
        etaSpanNext = findKnotSpan(xiEtaNext(2,1), Eta_surface, length(CP_surface(1,:,1)));
        
        %% 1vi.3. Compute the surface's basis functions at the left and the right segment's corners
        RPrevious = computeIGABasisFunctionsAndDerivativesForSurface...
            (xiSpanPrevious,p_surface,xiEtaPrevious(1,1),Xi_surface,etaSpanPrevious,q_surface,xiEtaPrevious(2,1),...
        Eta_surface,CP_surface,isNURBS_surface,0);
        RNext = computeIGABasisFunctionsAndDerivativesForSurface...
            (xiSpanNext,p_surface,xiEtaNext(1,1),Xi_surface,etaSpanNext,q_surface,xiEtaNext(2,1),...
            Eta_surface,CP_surface,isNURBS_surface,0);
        
        %% 1vi.4. Compute the Cartesian coordinates of the left and the right segment's corners
        SPrevious = computeCartesianCoordinatesOfAPointOnBSplineSurface...
            (xiSpanPrevious,p_surface,xiEtaPrevious(1,1),Xi_surface,etaSpanPrevious,q_surface,xiEtaPrevious(2,1),Eta_surface,CP_surface,RPrevious);
        SNext = computeCartesianCoordinatesOfAPointOnBSplineSurface...
            (xiSpanNext,p_surface,xiEtaNext(1,1),Xi_surface,etaSpanNext,q_surface,xiEtaNext(2,1),Eta_surface,CP_surface,RNext);
        
        %% 1vi.5. Plot the line segment
        plot3([SPrevious(1,1) SNext(1,1)],[SPrevious(2,1) SNext(2,1)],[SPrevious(3,1) SNext(3,1)],'--magenta');

        %% 1vi.6. Assign the right segment's corner to the left one
        xiEtaPrevious = xiEtaNext;
        
        %% 1vi.7. Update the line parameter
        lambda = lambda + dlambda;
    end
end

end
