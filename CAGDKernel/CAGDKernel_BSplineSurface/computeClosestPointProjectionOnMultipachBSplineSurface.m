function patchIDXiEta = ...
    computeClosestPointProjectionOnMultipachBSplineSurface...
    (BSplinePatches, P, epsProj, noXi, noEta, propNewtonRaphson,outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the closest point projection of a Point to a multipatch geometry.
% The point might have images into more than one patches depending on
% whether it is found on a common interface.
%
%              Input :
%     BSplinePatches : A cell array of all patches in the multipatch 
%                      geometry each of which containing:
%                        .p,.q : The polymomial orders of the patch along 
%                                xi and eta parametric directions
%                     .Xi,.Eta : The knot vectors of the patch along xi and
%                                eta parametric directions
%                          .CP : The set of coordinates and weights of the
%                                Control Points
%                     .isNURBS : Flag on whether the basis is a NURBS or a
%                                B-Spline
%                  P : Cartesian coordinates of the point to be projected
%            epsProj : Tolerance for the projection of the point on the
%                      multipatch geometry
%         noXi,noEta : Number of points used as initial guesses for the
%                      Newton-Raphson method
%  propNewtonRaphson : Properties of the Newton-Rapshon method for
%                      projection of the point on the NURBS multipatch
%                      geometry :
%                          .eps : Convergence tolerance
%                        .maxIt : Maximum number of Newton-Rapshon
%                                 iterations
%             outMsg : Allows printing onto command window if chosen as 
%                      'outputEnabled'
%
%       Output :
% patchIDXiEta : The id of the patch where the point was projected onto and
%                the corresponding parametric coordinates on the patch and
%                the Cartesian coordinates of the point at the parametric
%                location on the patch, 
%                patchIDXiEta = [patchID xi eta Px Py Pz]
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all patches
% ->
%    1i. Get the patch properties
%
%   1ii. Get number of points where to compute intial guesses for the projection of the nodes on the patch
%
%  1iii. Loop over all the initial guesses for the projection of the point
%
%   1iv. Update the output array counter if projection was found
%
%    1v. Fill up the output arrays accordingly
% <-
%
% 2. Check if any projection was found
%
% 3. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('___________________________________________________________________________\n');
    fprintf('###########################################################################\n');
    fprintf('Projection of a point onto a multipatch B-Spline surface has been initiated\n\n');
    fprintf('Projecting point P = (%d,%d,%d) \n',P(1,1),P(2,1),P(3,1));
    fprintf('Number of patches equals %d\n',length(BSplinePatches))
    fprintf('Projection tolerance equals %d \n',epsProj);
    fprintf('___________________________________________________________________________\n\n');

    % start measuring computational time
    tic;
end

%% 0. Read input

% Number of patches
noPatches = length(BSplinePatches);

% Initialize output arrays
patchIDXiEta = [];

% Initialize counter
counter = 0;

%% 1. Loop over all patches
for iPatches = 1:noPatches
    %% 1i. Get the patch properties
    p = BSplinePatches{iPatches}.p;
    Xi = BSplinePatches{iPatches}.Xi;
    q = BSplinePatches{iPatches}.q;
    Eta = BSplinePatches{iPatches}.Eta;
    CP = BSplinePatches{iPatches}.CP;
    isNURBS = BSplinePatches{iPatches}.isNURBS;
    isProjected = false;
    xi0 = .5; %Xi(1);
    eta0 = .5; %Eta(1);
    
    %% 1ii. Get number of points where to compute intial guesses for the projection of the nodes on the patch
    dxi = (Xi(length(Xi)) - Xi(1))/noXi;
    deta = (Eta(length(Eta)) - Eta(1))/noEta;
    
    %% 1iii. Loop over all the initial guesses for the projection of the point
    while xi0 < Xi(length(Xi)) && ~isProjected
        while eta0 < Eta(length(Eta)) && ~isProjected
            [xi,eta,PProj,isProjected,~] = ...
                computeNearestPointProjectionOnBSplineSurface...
                (P,p,Xi,q,Eta,CP,isNURBS,xi0,eta0,...
                propNewtonRaphson);
            if isProjected && norm(P - PProj) > epsProj
                isProjected = false;
            end
            eta0 = eta0 + deta;
        end
        xi0 = xi0 + dxi;
        eta0 = Eta(1);
    end
    
    %% 1iv. Update the output array counter if projection was found
    if isProjected
        counter = counter + 1;
        if strcmp(outMsg,'outputEnabled')
            fprintf('\t>> Point projected on patch %d with parametric image (xi,eta) = (%d,%d)\n',iPatches,xi,eta);
        end
    else
       continue; 
    end
    
    %% 1v. Fill up the output arrays accordingly
    patchIDXiEta(counter,:) = [iPatches xi eta PProj'];
end

%% 2. Check if any projection was found
if isempty(patchIDXiEta)
    warning('The point was not projected to any of the patches');
end

%% 3. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    fprintf('\nProjection of point onto multipatch B-Spline surface took %.2d seconds \n\n',computationalTime);
    fprintf('__________________Plotting Current Configuration Ended_____________________\n');
    fprintf('###########################################################################\n\n\n');
end

end
