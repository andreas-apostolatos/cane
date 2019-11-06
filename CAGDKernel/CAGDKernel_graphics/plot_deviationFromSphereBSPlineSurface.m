%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _______________________________________________________               %
%   _______________________________________________________               %
%                                                                         %
%   Technische Universität München                                        %
%   Lehrstuhl für Statik, Prof. Dr.-Ing. Kai-Uwe Bletzinger               %
%   _______________________________________________________               %
%   _______________________________________________________               %
%                                                                         %
%                                                                         %
%   Authors                                                               %
%   _______________________________________________________________       %
%                                                                         %
%   Dipl.-Math. Andreas Apostolatos    (andreas.apostolatos@tum.de)       %
%   Dr.-Ing. Roland Wüchner            (wuechner@tum.de)                  %
%   Prof. Dr.-Ing. Kai-Uwe Bletzinger  (kub@tum.de)                       %
%   _______________________________________________________________       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [relErrorInL2Domain,relErrorInL2WeakDBC,errorInL2Interface,index] = ...
    plot_deviationFromSphereBSPlineSurface...
    (BSplinePatches,dHat,connections,propSphere,int,graph,outMsg)
%% Function documentation
%
% Returns the relative error in the L2-norm of the difference of the given 
% B-Spline surface from an exact sphere for each of the patches as well as 
% the index of the current graph.
%
%                 Input :
%        BSplinePatches : Array of B-Spline patches each of which contains
%                             .p,.q : The polynomial orders of the B-Spline 
%                                     patch
%                          .Xi,.Eta : The knot vectors of the B-Spline 
%                                     patch
%                               .CP : The set of Control point coordinates 
%                                     and weights of the B-Spline patch
%                       .EFTPatches : The freedom tables for each B-Spline 
%                                     patch in the multipatch geometry
%           connections : Define the connection between the patches:
%                            .No : Number of connections
%                     .xiEtaCoup : [patchID1 patchID2 xi12 eta12 xi21 eta21
%                                   ...      ...    ...   ...  ...   ...]
%                  dHat : The displacement field corresponding to the 
%                         multipatch geometry, if chosen as a string the 
%                         variable CPd of the B-Spline patch is chosen
%            propSphere : Analytical data for the analytical description of 
%                         the sphere :
%                             .center : The center of the exact sphere
%                             .radius : The radius of the exact sphere
%                   int : Quadrature scheme for the computation of the L2-
%                         norm
%                             .type : 'default' or 'user'
%                           .noGPXi : Number of Gauss points for the 
%                                     integration in the xi-direction if 
%                                     int.type = 'user'
%                          .noGPEta : Number of Gauss points for the 
%                                     integration in the eta-direction if 
%                                     int.type = 'user'
%                 graph : On the graphics. If graph is chosen as 
%                         'undefined' no graphics are produced
%                outMsg : Enables outputting information onto the command 
%                         window when chosen as 'outputEnabled'
%
%                Output :
%    relErrorInL2Domain : The relative error in the L2 norm of the domain
%                         for each patch [noPatches,1]
%   relErrorInL2WeakDBC : The relative error in the L2-norm from the 
%                         satisfaction of the weakly applied Dirichlet
%                         boundary condition [noWeakDBC,IDPatch,IDCnd]
%    errorInL2Interface : The relative error in the L2-norm for the jump
%                         accross the interface between the neighbouring
%                         patches [noConnections,1]
%                 index : The index of the current graph
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all patches in the B-Spline patch geometry
% ->
%    1i. Get the patch parameters
%
%   1ii. Compute the displaced Control Points of the patch
%
%  1iii. Get the Gauss points and weights of the quadrature for the computation of the L2-norm
%
%   1iv. Initialize the array of the deviation from a sphere for the plot
%
%    1v. Loop over all the elements
%    ->
%        1v.1. Loop over all Gauss points
%
%        1v.2. Loop over all Gauss points
%        ->
%              1v.2i. Compute the NURBS coordinates xi,eta of the Gauss Point coordinates in the bi-unit interval [-1, 1]
%
%             1v.2ii. Compute the NURBS basis functions and their first derivatives at the Gauss Point
%
%            1v.2iii. Compute the covariant base vectors of the reference configuration at the Gauss point
%
%             1v.2iv. Compute the surface normal of the reference configuration at the Gauss point (third covariant base vector not normalized)
%
%              1v.2v. Compute the legth of G3Tilde at the Gauss point (= area dA of the undeformed configuration)
%
%             1v.2vi. Compute the point on the deformed B-Spline surface
%
%            1v.2vii. Compute the direction vector of the computed point on the B-Spline surface
%
%           1v.2viii. Find the corresponding spherical coordinates with respect to the computed direction vector
%
%             1v.2ix. Compute the coordinates of the point on the sphere
%
%              1v.2x. Compute the element area on the Gauss Point
%
%             1v.2xi. Compute the L2-norm of the difference from a sphere and the norm of the exact sphere
%        <-
%    <-
%
%   1vi. Loop over the boundaries of the patch where weakly Dirichlet boundary conditions are applied
%   ->
%        1vi.1. Assign the numbering of the corresponding patch and the corresponding Dirichlet boundary condition
%
%        1vi.2. Get the boundary of the patch where weakly Dirichlet boundary conditions are applied
%
%        1vi.3. Check along which parametric direction the Dirichlet boundary conditions are weakly applied
%
%        1vi.4. Choose an boundary integration rule
%
%        1vi.5. Loop over the elements on the coupling interface
%        ->
%               1vi.5i. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
%
%              1vi.5ii. Loop over the Gauss points
%              ->
%                       1vi.5ii.1. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
%
%                       1vi.5ii.2. Compute the NURBS basis functions
%
%                       1vi.5ii.3. Compute the Cartesian coordinates of the point on the boundary
%
%                       1vi.5ii.4. Compute the expected value for the form-found shape
%
%                       1vi.5ii.5. Compute the covariant base vectors of the reference configuration
%
%                       1vi.5ii.6. Compute the determinant of the Jacobian of the transformation from the physical to the parent domain on the GP
%
%                       1vi.5ii.7. Compute the element length at the GP
%
%                       1vi.5ii.8. Compute the difference of the displacements from each patch and the the norm of the expected solution at the Gauss point and add the contribution
%               <-
%
%             1vi.5iii. Update the counter of the Dirichlet conditions which are weakly applied
%        <-
%   <-
%
%  1vii. Initialize the eta parametric coordinate and the counter in eta-direction
%
% 1viii. Loop over all the sampling points in -eta direction
% ->
%        1viii.1. Find the span in the eta-direction
%
%        1viii.2. Initialize the xi parametric coordinate and the counter in xi-direction
%
%        1viii.3. Loop over all the sampling points in -xi direction
%        ->
%                 1viii.3i. Find the span in xi-direction
%
%                1viii.3ii. Compute the IGA basis functions and their derivatives
%
%               1viii.3iii. Compute the point on the deformed B-Spline surface
%
%                1viii.3iv. Compute the direction vector of the computed point on the B-Spline surface
%
%                 1viii.3v. Find the corresponding spherical coordinates with respect to the computed direction vector
%
%                1viii.3vi. Compute the coordinates of the point on the sphere
%
%               1viii.3vii. Compute the deviation from the exact sphere at the samping point
%
%              1viii.3viii. Update the parametric coordinate and the counter in xi-direction
%        <-
%
%        1viii.4. Update the parametric coordinate in eta-direction and the counter
% <-
%
%   1ix. Plot the deformed surface together with the difference from an exact sphere
% <-
%
% 2. Loop over all the connections in the multipatch geometry
% ->
%    2i. Get the patch IDs
%
%   2ii. Get the coupling boundaries for both patches
%
%  2iii. Determine the interface orientation
%
%   2iv. Recover the B-Spline patch parameters
%
%    2v. Get the displacement vector for each patch and compute the displaced Control Points
%
%   2vi. Get the running and the fixed parameters on the patch interface and the coupling region
%
%  2vii. Compute the merged knot vector from both patches over the interface
%
% 2viii. Choose an interface integration rule
%
%   2ix. Loop over the elements on the coupling interface
%   ->
%        2ix.1. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
%
%        2ix.2. Loop over the Gauss points
%        ->
%               2ix.2i. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
%
%              2ix.2ii. Compute the NURBS basis functions
%
%             2ix.2iii. Compute the Cartesian components of the points on the interface boundary
%
%              2ix.2iv. Compute the determinant of the Jacobian of the transformation from the physical to the parent domain on the GP
%
%               2ix.2v. Compute the element length at the GP
%
%              2ix.2vi. Compute the difference of the displacements from each neighbouring patch
%        <-
%   <-
% <-
%
% 3. Assign graphic properties
%
% 4. Compute the relative error in the L2-norm for the deviation from the exact sphere
%
% 5. Appendix
% 
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('_________________________________________________________________________\n');
    fprintf('#########################################################################\n');
    fprintf('Plotting of the distribution of the deviation from an exact sphere over a\n');
    fprintf('multipatch B-Spline surface and computation of the corresponding L2-norms\n');
    fprintf('has been initiated \n');
    fprintf('_________________________________________________________________________\n\n');
    
    % start measuring computational time
    tic;
end
%% 0. Read input

% Number of points for the visualization of the deviation from the sphere
xiGrid = 49;
etaGrid = 49;

% Assign a tolerance value
tol = 1e-9;

% Number of patches
noPatches = length(BSplinePatches);

% Check the number of the Dirichlet boundaries where weakly imposed
% conditions are assumed
noCnd = 0;
for iPatches = 1:noPatches
    if isfield(BSplinePatches{iPatches},'weakDBC')
        if ~ischar(BSplinePatches{iPatches}.weakDBC)
            if isfield(BSplinePatches{iPatches}.weakDBC,'noCnd')
                if ~ischar(BSplinePatches{iPatches}.weakDBC.noCnd)
                    for iCnd = 1:BSplinePatches{iPatches}.weakDBC.noCnd
                        noCnd = noCnd + 1;
                    end
                end
            end
        end
    end
end

% Initialize error-related arrays for the domain error
errorInL2Domain = zeros(noPatches,1);
exactInL2Domain = zeros(noPatches,1);
relErrorInL2Domain = zeros(noPatches,1);

% Initialize error-related arrays for the error along the Dirichlet
% boundaries
errorInL2WeakDBC = zeros(noCnd,3);
exactInL2WeakDBC = zeros(noCnd,3);
relErrorInL2WeakDBC = zeros(noCnd,3);

% Initialize the error-related arrays for the error along the interface
% boundaries
if ~isempty(connections)
    errorInL2Interface = zeros(connections.No,3);
else
    errorInL2Interface = 'undefined';
end

% Initialize the counter of the boundaries where weakly Dirichlet boundary conditions are applied
counterCnd = 1;

% Initialize handle to the figure
if ~strcmp(graph,'undefined')
    figure(graph.index)
else
    index = 'undefined';
end

%% 1. Loop over all patches in the B-Spline patch geometry
for iPatches = 1:noPatches
    %% 1i. Get the patch parameters
    p = BSplinePatches{iPatches}.p;
    q = BSplinePatches{iPatches}.q;
    Xi = BSplinePatches{iPatches}.Xi;
    Eta = BSplinePatches{iPatches}.Eta;
    CP = BSplinePatches{iPatches}.CP;
    isNURBS = BSplinePatches{iPatches}.isNURBS;
    weakDBC = BSplinePatches{iPatches}.weakDBC;
    nxi = length(CP(:,1,1));
    neta = length(CP(1,:,1));
    mxi = length(Xi);
    mxiUnique = length(unique(Xi));
    meta = length(Eta);
    metaUnique = length(unique(Eta));
    if ~ischar(dHat)
        EFTPatches = BSplinePatches{iPatches}.EFTPatches;
    end
    dxi = (Xi(mxi)-Xi(1))/(etaGrid-1);
    deta = (Eta(meta)-Eta(1))/(xiGrid-1);
    
    %% 1ii. Compute the displaced Control Points of the patch
    if ~ischar(dHat)
        CPd = computeDisplacedControlPointsForIGAKirchhoffLoveShell...
            (CP,dHat(EFTPatches));
    else
        CPd = BSplinePatches{iPatches}.CPd;
    end
    
    %% 1iii. Get the Gauss points and weights of the quadrature for the computation of the L2-norm
    if strcmp(int.type,'default')
        noGPXi = ceil((p + 1)/2);
        noGPEta = ceil((q + 1)/2);
        isQuadratureUser = false;
    elseif strcmp(int.type,'user')
        noGPXi = int.noGPXi;
        noGPEta = int.noGPEta;
        [GPXi,GWXi] = getGaussPointsAndWeightsOverUnitDomain(noGPXi);
        [GPEta,GWEta] = getGaussPointsAndWeightsOverUnitDomain(noGPEta);
        isQuadratureUser = true;
    else
        error('Define the type of quadrature rule');
    end
    if ~isQuadratureUser
        [GPXi,GWXi] = getGaussPointsAndWeightsOverUnitDomain(noGPXi);
        [GPEta,GWEta] = getGaussPointsAndWeightsOverUnitDomain(noGPEta);
    end
    
	%% 1iv. Initialize the array of the deviation from a sphere for the plot
    clear deviationFromSphere X;
    deviationFromSphere = zeros(noGPXi*mxiUnique,noGPEta*metaUnique);
    P = zeros(noGPXi*mxiUnique,noGPEta*metaUnique,3);
    
    %% 1v. Loop over all the elements
    for j = q+1:meta-q-1
        if Eta(j+1) ~= Eta(j)
            for i = p+1:mxi-p-1
                if Xi(i+1) ~= Xi(i)
                    %% 1v.1. Loop over all Gauss points
                    %
                    %         | xi_i+1 - xi_i                    |
                    %         | -------------            0       |
                    %         |        2                         |
                    %  xi,u = |                                  |
                    %         |                  eta_j+1 - eta_j |
                    %         |        0         --------------- |
                    %         |                          2       |
                    detJxiu = (Xi(i+1)-Xi(i))*(Eta(j+1)-Eta(j))/4;

                    %% 1v.2. Loop over all Gauss points
                    for iGPEta = 1:noGPEta
                        for iGPXi = 1:noGPXi
                            %% 1v.2i. Compute the NURBS coordinates xi,eta of the Gauss Point coordinates in the bi-unit interval [-1, 1]
                            xi = ( Xi(i+1)+Xi(i) + GPXi(iGPXi)*(Xi(i+1)-Xi(i)) )/2;
                            eta = ( Eta(j+1)+Eta(j) + GPEta(iGPEta)*(Eta(j+1)-Eta(j)) )/2;

                            %% 1v.2ii. Compute the NURBS basis functions and their first derivatives at the Gauss Point
                            dR = computeIGABasisFunctionsAndDerivativesForSurface...
                                (i,p,xi,Xi,j,q,eta,Eta,CP,isNURBS,1);

                            %% 1v.2iii. Compute the covariant base vectors of the reference configuration at the Gauss point
                            [A1,A2] = computeBaseVectorsAndDerivativesForBSplineSurface...
                                (i,p,j,q,CP,0,dR);

                            %% 1v.2iv. Compute the surface normal of the reference configuration at the Gauss point (third covariant base vector not normalized)
                            A3Tilde = cross(A1(:,1),A2(:,1));

                            %% 1v.2v. Compute the legth of G3Tilde at the Gauss point (= area dA of the undeformed configuration)
                            dA = norm(A3Tilde);

                            %% 1v.2vi. Compute the point on the deformed B-Spline surface
                            X = computeCartesianCoordinatesOfAPointOnBSplineSurface...
                                (i,p,xi,Xi,j,q,eta,Eta,CPd,dR(:,1));

                            %% 1v.2vii. Compute the direction vector of the computed point on the B-Spline surface
                            directionVct = X - propSphere.center;

                            %% 1v.2viii. Find the corresponding spherical coordinates with respect to the computed direction vector
                            theta = atan(directionVct(2,1)/directionVct(1,1));
                            phi = acos(directionVct(3,1)/norm(directionVct));

                            %% 1v.2ix. Compute the coordinates of the point on the sphere
                            XExact = [propSphere.radius*cos(theta)*sin(phi)
                                       propSphere.radius*sin(theta)*sin(phi)
                                       propSphere.radius*cos(phi)];

                            %% 1v.2x. Compute the element area on the Gauss Point
                            elementAreaOnGP = dA*detJxiu*GWXi(iGPXi)*GWEta(iGPEta);

                            %% 1v.2xi. Compute the L2-norm of the difference from a sphere and the norm of the exact sphere
                            errorInL2Domain(iPatches,1) = errorInL2Domain(iPatches,1) + ...
                                norm(abs(X) - abs(XExact))^2*elementAreaOnGP;
                            exactInL2Domain(iPatches,1) = exactInL2Domain(iPatches,1) + ...
                                norm(XExact)^2*elementAreaOnGP;
                        end
                    end
                end
            end
        end
    end
    
    %% 1vi. Loop over the boundaries of the patch where weakly Dirichlet boundary conditions are applied
    if ~ischar(BSplinePatches{iPatches}.weakDBC)
        if isfield(BSplinePatches{iPatches}.weakDBC,'noCnd')
            if ~ischar(BSplinePatches{iPatches}.weakDBC.noCnd)
                for iCnd = 1:BSplinePatches{iPatches}.weakDBC.noCnd
                    %% 1vi.1. Assign the numbering of the corresponding patch and the corresponding Dirichlet boundary condition
                    errorInL2WeakDBC(counterCnd,2) = iPatches;
                    errorInL2WeakDBC(counterCnd,3) = iCnd;

                    %% 1vi.2. Get the boundary of the patch where weakly Dirichlet boundary conditions are applied
                    xiWeakDBC = weakDBC.xiExtension{iCnd};
                    etaWeakDBC = weakDBC.etaExtension{iCnd};

                    %% 1vi.3. Check along which parametric direction the Dirichlet boundary conditions are weakly applied
                    if etaWeakDBC(1) == etaWeakDBC(2)
                        % Coupled region in xi-direction
                        weakDBCRegion = xiWeakDBC;

                        % Find the correct spans for the coupled region
                        spanStart = findKnotSpan(weakDBCRegion(1),Xi,nxi);
                        spanEnd = findKnotSpan(weakDBCRegion(2),Xi,nxi) + 1;

                        % Corresponding to the coupled region knot span
                        weakDBCRegionOnKnotVector = Xi(spanStart:spanEnd);

                        % Fixed parameter on the parametric net
                        eta = etaWeakDBC(1);

                        % Find the span where xiEta it lies in
                        etaSpan = findKnotSpan(eta,Eta,neta);

                        % Flag on whether the coupling line is over xi
                        isOnXi = true;
                    else
                        % Coupled region in eta-direction
                        weakDBCRegion = etaWeakDBC;

                        % Find the correct spans for the coupled region
                        spanStart = findKnotSpan(weakDBCRegion(1),Eta,neta);   
                        spanEnd = findKnotSpan(weakDBCRegion(2),Eta,neta) + 1;   

                        % Corresponding to the coupled region knot span
                        weakDBCRegionOnKnotVector = Eta(spanStart:spanEnd);

                        % Fixed parameter on the parametric net
                        xi = xiWeakDBC(1);

                        % Find the span where uv it lies in
                        xiSpan = findKnotSpan(xi,Xi,nxi);

                        % Flag on whether the coupling line is over eta
                        isOnXi = false;
                    end

                    %% 1vi.4. Choose an boundary integration rule
                    if strcmp(int.type,'default')
                        if isOnXi
                            noGPs = p + 1;
                        else
                            noGPs = q + 1;
                        end
                    elseif strcmp(int.type,'user')
                        noGPs = int.noGPs;
                    end
                    [GP,GW] = getGaussPointsAndWeightsOverUnitDomain(noGPs);

                    %% 1vi.5. Loop over the elements on the coupling interface
                    for iElmts = 1:length(weakDBCRegionOnKnotVector)-1
                        if weakDBCRegionOnKnotVector(iElmts) ~= weakDBCRegionOnKnotVector(iElmts + 1)
                            %% 1vi.5i. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
                            detJxizeta = (weakDBCRegionOnKnotVector(iElmts + 1) - weakDBCRegionOnKnotVector(iElmts))/2;

                            %% 1vi.5ii. Loop over the Gauss points
                            for iGPs = 1:noGPs
                                %% 1vi.5ii.1. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
                                xiEta = ((1-GP(iGPs))*weakDBCRegionOnKnotVector(iElmts) + ...
                                    (1+GP(iGPs))*weakDBCRegionOnKnotVector(iElmts + 1))/2;

                                %% 1vi.5ii.2. Compute the NURBS basis functions
                                if isOnXi
                                    xi = xiEta;
                                    xiSpan = findKnotSpan(xi,Xi,nxi);
                                else
                                    eta = xiEta;
                                    etaSpan = findKnotSpan(eta,Eta,neta);
                                end
                                dR = computeIGABasisFunctionsAndDerivativesForSurface...
                                    (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,1);

                                %% 1vi.5ii.3. Compute the Cartesian coordinates of the point on the boundary
                                X = computeCartesianCoordinatesOfAPointOnBSplineSurface...
                                    (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CPd,dR(:,1));

                                %% 1vi.5ii.4. Compute the expected value for the form-found shape

                                % Compute the direction vector
                                directionVct = X - propSphere.center;

                                % Find the corresponding spherical coordinates with respect to the computed direction vector
                                theta = atan(directionVct(2,1)/directionVct(1,1));
                                phi = acos(directionVct(3,1)/norm(directionVct));

                                % Compute the coordinates of the point on the sphere
                                XExact = [propSphere.radius*cos(theta)*sin(phi)
                                          propSphere.radius*sin(theta)*sin(phi)
                                          propSphere.radius*cos(phi)];

                                %% 1vi.5ii.5. Compute the covariant base vectors of the reference configuration
                                [A1,A2] = computeBaseVectorsAndDerivativesForBSplineSurface...
                                    (xiSpan,p,etaSpan,q,CP,0,dR);

                                %% 1vi.5ii.6. Compute the determinant of the Jacobian of the transformation from the physical to the parent domain on the GP
                                if isOnXi
                                    detJxxi = norm(A1(:,1));
                                else
                                    detJxxi = norm(A2(:,1));
                                end

                                %% 1vi.5ii.7. Compute the element length at the GP
                                elementLengthOnGP = detJxxi*detJxizeta*GW(iGPs);

                                %% 1vi.5ii.8. Compute the difference of the displacements from each patch and the the norm of the expected solution at the Gauss point and add the contribution
                                errorInL2WeakDBC(counterCnd,1) = ...
                                    errorInL2WeakDBC(counterCnd,1) + norm(abs(X) - abs(XExact))*elementLengthOnGP;
                                exactInL2WeakDBC(counterCnd,1) = ...
                                    exactInL2WeakDBC(counterCnd,1) + norm(abs(XExact) - abs(XExact))*elementLengthOnGP;
                            end
                        end
                    end

                    %% 1vi.5iii. Update the counter of the Dirichlet conditions which are weakly applied
                    counterCnd = counterCnd + 1;
                end
            end
        end
    end
    
    if ~strcmp(graph,'undefined')
        %% 1vii. Initialize the eta parametric coordinate and the counter in eta-direction
        eta = Eta(1);
        etaCounter = 1;  

        %% 1viii. Loop over all the sampling points in -eta direction
        while eta <= Eta(meta)+tol
            %% 1viii.1. Find the span in the eta-direction
            etaSpan = findKnotSpan(eta,Eta,neta);

            %% 1viii.2. Initialize the xi parametric coordinate and the counter in xi-direction
            xi = Xi(1);
            xiCounter = 1;

            %% 1viii.3. Loop over all the sampling points in -xi direction
            while xi <= Xi(mxi)+tol
                %% 1viii.3i. Find the span in xi-direction
                xiSpan = findKnotSpan(xi,Xi,nxi);

                %% 1viii.3ii. Compute the IGA basis functions and their derivatives
                R = computeIGABasisFunctionsAndDerivativesForSurface...
                    (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,0);

                %% 1viii.3iii. Compute the point on the deformed B-Spline surface
                P(xiCounter,etaCounter,1:3) = ...
                    computeCartesianCoordinatesOfAPointOnBSplineSurface...
                    (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CPd,R);

                %% 1viii.3iv. Compute the direction vector of the computed point on the B-Spline surface
                directionVct = abs(squeeze(P(xiCounter,etaCounter,:)) - propSphere.center);

                %% 1viii.3v. Find the corresponding spherical coordinates with respect to the computed direction vector
                theta = atan(directionVct(2,1)/directionVct(1,1));
                phi = acos(directionVct(3,1)/norm(directionVct));

                %% 1viii.3vi. Compute the coordinates of the point on the sphere
                XExact = [propSphere.radius*cos(theta)*sin(phi)
                          propSphere.radius*sin(theta)*sin(phi)
                          propSphere.radius*cos(phi)];

                %% 1viii.3vii. Compute the deviation from the exact sphere at the samping point
                unitVct = (XExact - propSphere.center)/norm(XExact - propSphere.center);
                deviationFromSphere(xiCounter,etaCounter) = ...
                    (abs(squeeze(P(xiCounter,etaCounter,:))) - abs(XExact))'*unitVct;

                %% 1viii.3viii. Update the parametric coordinate and the counter in xi-direction
                xi = xi + dxi;
                xiCounter = xiCounter + 1;
             end

            %% 1viii.4. pdate the parametric coordinate in eta-direction and the counter
            eta = eta + deta;
            etaCounter = etaCounter + 1;
        end

        %% 1ix. Plot the deformed surface together with the difference from an exact sphere
        surf(P(:,:,1),P(:,:,2),P(:,:,3),deviationFromSphere(:,:));
        hold on;
        plot_knotsForBSplineSurfaceOnCartesianSpace(p,q,Xi,Eta,CPd,isNURBS,false,xiGrid,etaGrid);
    end
end
if ~strcmp(graph,'undefined')
    hold off;
end

%% 2. Loop over all the connections in the multipatch geometry
if ~isempty(connections)
    for iConnections = 1:connections.No
        %% 2i. Get the patch IDs

        % Patch I :
        % _________

        IDI = connections.xiEtaCoup(iConnections,1);
        errorInL2Interface(iConnections,2) = IDI;
        errorInL2Interface(iConnections,2) = IDI;

        % Patch J :
        % _________

        IDJ = connections.xiEtaCoup(iConnections,2);
        errorInL2Interface(iConnections,3) = IDJ;
        errorInL2Interface(iConnections,3) = IDJ;

        %% 2ii. Get the coupling boundaries for both patches

        % Patch I :
        % _________

        BSplinePatches{IDI}.xicoup = connections.xiEtaCoup(iConnections,3:4);
        BSplinePatches{IDI}.etacoup = connections.xiEtaCoup(iConnections,5:6);

        % Patch J :
        % _________

        BSplinePatches{IDJ}.xicoup = connections.xiEtaCoup(iConnections,7:8);
        BSplinePatches{IDJ}.etacoup = connections.xiEtaCoup(iConnections,9:10);

        %% 2iii. Determine the interface orientation
            haveSameOrientation = findSubdomainInterfaceOrientation...
                (BSplinePatches{IDI}.p,BSplinePatches{IDI}.Xi,BSplinePatches{IDI}.q,BSplinePatches{IDI}.Eta,BSplinePatches{IDI}.CP,BSplinePatches{IDI}.isNURBS,BSplinePatches{IDI}.xicoup,BSplinePatches{IDI}.etacoup,...
                BSplinePatches{IDJ}.p,BSplinePatches{IDJ}.Xi,BSplinePatches{IDJ}.q,BSplinePatches{IDJ}.Eta,BSplinePatches{IDJ}.CP,BSplinePatches{IDJ}.isNURBS,BSplinePatches{IDJ}.xicoup,BSplinePatches{IDJ}.etacoup);

        %% 2iv. Recover the B-Spline patch parameters

        % Patch I :
        % _________

        pI = BSplinePatches{IDI}.p;
        qI = BSplinePatches{IDI}.q;
        XiI = BSplinePatches{IDI}.Xi;
        EtaI = BSplinePatches{IDI}.Eta;
        CPI = BSplinePatches{IDI}.CP;
        isNURBSI = BSplinePatches{IDI}.isNURBS;
        xicoupI = BSplinePatches{IDI}.xicoup;
        etacoupI = BSplinePatches{IDI}.etacoup;
        nxiI = length(CPI(:,1,1));
        netaI = length(CPI(1,:,1));

        % Patch J :
        % _________

        pJ = BSplinePatches{IDJ}.p;
        qJ = BSplinePatches{IDJ}.q;
        XiJ = BSplinePatches{IDJ}.Xi;
        EtaJ = BSplinePatches{IDJ}.Eta;
        CPJ = BSplinePatches{IDJ}.CP;
        isNURBSJ = BSplinePatches{IDJ}.isNURBS;
        xicoupJ = BSplinePatches{IDJ}.xicoup;
        etacoupJ = BSplinePatches{IDJ}.etacoup;
        nxiJ = length(CPJ(:,1,1));
        netaJ = length(CPJ(1,:,1));

        %% 2v. Get the displacement vector for each patch and compute the displaced Control Points

        % For patch I :
        % _____________

        if ~ischar(dHat)
            CPdI = computeDisplacedControlPointsForIGAKirchhoffLoveShell...
                (CPI,dHatI);
        else
            CPdI = BSplinePatches{IDI}.CPd;
        end

        % For patch J :
        % _____________

        if ~ischar(dHat)
            CPdJ = computeDisplacedControlPointsForIGAKirchhoffLoveShell...
                (CPJ,dHatJ);
        else
            CPdJ = BSplinePatches{IDJ}.CPd;
        end

        %% 2vi. Get the running and the fixed parameters on the patch interface and the coupling region

        % For patch I :
        % _____________

        if etacoupI(1) == etacoupI(2)
            % Coupled region in xi-direction
            couplingRegionI = xicoupI;

            % Find the correct spans for the coupled region
            spanStartI = findKnotSpan(couplingRegionI(1),XiI,nxiI);
            spanEndI = findKnotSpan(couplingRegionI(2),XiI,nxiI)+1;

            % Corresponding to the coupled region knot span
            couplingRegionOnKnotVectorI = XiI(spanStartI:spanEndI);

            % Fixed parameter on the parametric net
            etaI = etacoupI(1);

            % Find the span where xiEta it lies in
            etaSpanI = findKnotSpan(etaI,EtaI,netaI);

            % Flag on whether the coupling line is over xi
            isOnXiI = true;
        else
            % Coupled region in eta-direction
            couplingRegionI = etacoupI;

            % Find the correct spans for the coupled region
            spanStartI = findKnotSpan(couplingRegionI(1),EtaI,netaI);   
            spanEndI = findKnotSpan(couplingRegionI(2),EtaI,netaI)+1;   

            % Corresponding to the coupled region knot span
            couplingRegionOnKnotVectorI = EtaI(spanStartI:spanEndI);

            % Fixed parameter on the parametric net
            xiI = xicoupI(1);

            % Find the span where uv it lies in
            xiSpanI = findKnotSpan(xiI,XiI,nxiI);

            % Flag on whether the coupling line is over eta
            isOnXiI = false;
        end

        % For patch J :
        % _____________

        if etacoupJ(1) == etacoupJ(2)
            % Coupled region in xi-direction
            couplingRegionJ = xicoupJ;

            % Find the correct spans for the coupled region
            spanStartJ = findKnotSpan(couplingRegionJ(1),XiJ,nxiJ);   
            spanEndJ = findKnotSpan(couplingRegionJ(2),XiJ,nxiJ)+1; 

            % Corresponding to the coupled region knot span
            couplingRegionOnKnotVectorJ = XiJ(spanStartJ:spanEndJ);

            % Fixed parameter on the parametric net
            etaJ = etacoupJ(1);

            % Find the span where xiEta it lies in
            etaSpanJ = findKnotSpan(etaJ,EtaJ,netaJ);

            % Flag on whether the coupling line is over xi
            isOnXiJ = true;
        else
            % Coupled region in eta-direction
            couplingRegionJ = etacoupJ;

            % Find the correct spans for the coupled region
            spanStartJ = findKnotSpan(couplingRegionJ(1),EtaJ,netaJ);   
            spanEndJ = findKnotSpan(couplingRegionJ(2),EtaJ,netaJ)+1;

            % Corresponding to the coupled region knot span
            couplingRegionOnKnotVectorJ = EtaJ(spanStartJ:spanEndJ);

            % Fixed parameter on the parametric net
            xiJ = xicoupJ(1);

            % Find the span where uv it lies in
            xiSpanJ = findKnotSpan(xiJ,XiJ,nxiJ);

            % Flag on whether the coupling line is over eta
            isOnXiJ = false;
        end

        %% 2vii. Compute the merged knot vector from both patches over the interface

        % Merge the two knot vectors into one for integration purposes:
        couplingRegionOnKnotVector = mergesorted(couplingRegionOnKnotVectorI,couplingRegionOnKnotVectorJ);

        % Delete double entries
        couplingRegionOnKnotVector = unique(couplingRegionOnKnotVector);

        %% 2viii. Choose an interface integration rule
        if strcmp(int.type,'default')
            if isOnXiI
                pDegreeI = pI + 1;
            else
                pDegreeI = qI + 1;
            end
            if isOnXiJ
                pDegreeJ = pJ + 1;
            else
                pDegreeJ = qJ + 1;
            end
            noGPs = ceil((pDegreeI + pDegreeJ + 1)/2);
        elseif strcmp(int.type,'user')
            noGPs = int.nGP;
        end
        [GP,GW] = getGaussPointsAndWeightsOverUnitDomain(noGPs);

        %% 2ix. Loop over the elements on the coupling interface
        for iElmts = 1:length(couplingRegionOnKnotVector)-1
            if couplingRegionOnKnotVector(iElmts) ~= couplingRegionOnKnotVector(iElmts+1)
                %% 2ix.1. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
                detJxizeta = (couplingRegionOnKnotVector(iElmts+1) - couplingRegionOnKnotVector(iElmts))/2;

                %% 2ix.2. Loop over the Gauss points
                for iGPs = 1:noGPs
                    %% 2ix.2i. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
                    xiEta = ((1-GP(iGPs))*couplingRegionOnKnotVector(iElmts) + ...
                        (1+GP(iGPs))*couplingRegionOnKnotVector(iElmts+1))/2;

                    %% 2ix.2ii. Compute the NURBS basis functions

                    % For patch I :
                    % _____________

                    if isOnXiI
                        xiI = xiEta;
                        xiSpanI = findKnotSpan(xiI,XiI,nxiI);
                    else
                        etaI = xiEta;
                        etaSpanI = findKnotSpan(etaI,EtaI,netaI);
                    end
                    dRI = computeIGABasisFunctionsAndDerivativesForSurface...
                        (xiSpanI,pI,xiI,XiI,etaSpanI,qI,etaI,EtaI,CPI,isNURBSI,1);

                    % For patch J :
                    % _____________

                    if isOnXiJ
                        xiJ = xiEta;
                        if ~haveSameOrientation
                            xiJ = XiJ(length(XiJ)) - xiJ;
                        end
                        xiSpanJ = findKnotSpan(xiJ,XiJ,nxiJ);
                    else
                        etaJ = xiEta;
                        if ~haveSameOrientation
                            etaJ = EtaJ(length(EtaJ)) - etaJ;
                        end
                        etaSpanJ = findKnotSpan(etaJ,EtaJ,netaJ);
                    end
                    dRJ = computeIGABasisFunctionsAndDerivativesForSurface...
                        (xiSpanJ,pJ,xiJ,XiJ,etaSpanJ,qJ,etaJ,EtaJ,CPJ,isNURBSJ,1);

                    %% 2ix.2iii. Compute the Cartesian components of the points on the interface boundary

                    % For patch I :
                    % _____________

                    XI = computeCartesianCoordinatesOfAPointOnBSplineSurface...
                        (xiSpanI,pI,xiI,XiI,etaSpanI,qI,etaI,EtaI,CPdI,dRI(:,1));

                    % For patch J :
                    % _____________

                    XJ = computeCartesianCoordinatesOfAPointOnBSplineSurface...
                        (xiSpanJ,pJ,xiJ,XiJ,etaSpanJ,qJ,etaJ,EtaJ,CPdJ,dRJ(:,1));

                    %% 2ix.2iv. Compute the determinant of the Jacobian of the transformation from the physical to the parent domain on the GP
                    [A1I,A2I] = computeBaseVectorsAndDerivativesForBSplineSurface...
                        (xiSpanI,pI,etaSpanI,qI,CPI,0,dRI);
                    if isOnXiI
                        detJxxi = norm(A1I(:,1));
                    else
                        detJxxi = norm(A2I(:,1));
                    end

                    %% 2ix.2v. Compute the element length at the GP
                    elementLengthOnGP = detJxxi*detJxizeta*GW(iGPs);

                    %% 2ix.2vi. Compute the difference of the displacements from each neighbouring patch
                    errorInL2Interface(iConnections,1) = ...
                        errorInL2Interface(iConnections,1) + norm(XI - XJ)*elementLengthOnGP;
                end
            end
        end

    end
end

%% 3. Assign graphic properties
if ~strcmp(graph,'undefined')
    shading interp;
    title('Distribution of the deviation from an exact sphere');
    colormap('jet');
    colorbar;
    axis equal;
    xlabel('x','FontSize',14);
    ylabel('y','FontSize',14);
    zlabel('z','FontSize',14);
    index = graph.index + 1;
end

%% 4. Compute the relative error in the L2-norm for the deviation from the exact sphere
relErrorInL2WeakDBC(:,2:3) = errorInL2WeakDBC(:,2:3);
if ~isempty(connections)
    errorInL2Interface(:,2:3) = errorInL2Interface(:,2:3);
end
counterCnd = 1; 
for iPatches = 1:noPatches
    relErrorInL2Domain(iPatches,1) = ...
        sqrt(errorInL2Domain(iPatches,1))/sqrt(exactInL2Domain(iPatches,1));
    if ~ischar(BSplinePatches{iPatches}.weakDBC)
        if isfield(BSplinePatches{iPatches}.weakDBC,'noCnd')
            if ~ischar(BSplinePatches{iPatches}.weakDBC.noCnd)
                for iCnd = 1:BSplinePatches{iPatches}.weakDBC.noCnd
                    if norm(exactInL2WeakDBC(counterCnd,1)) ~= 0
                        relErrorInL2WeakDBC(counterCnd,1) = ...
                            sqrt(errorInL2WeakDBC(counterCnd,1))/sqrt(exactInL2WeakDBC(counterCnd,1));
                    else
                        relErrorInL2WeakDBC(counterCnd,1) = ...
                            sqrt(errorInL2WeakDBC(counterCnd,1));
                    end
                    counterCnd = counterCnd + 1;
                end
            end
        end
    end
end
if ~isempty(connections)
    for iConnections = 1:connections.No
        errorInL2Interface(iConnections,1) = sqrt(errorInL2Interface(iConnections,1));
    end
end

%% 5. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;
    if ~strcmp(graph,'undefined')
        fprintf('Plotting and computing the error took %.2d seconds \n\n',computationalTime);
        fprintf('__________Plotting And Computing The Error Distribution Ended____________\n');
        fprintf('#########################################################################\n\n\n');
    else
        fprintf('Plotting and computing the error took %.2d seconds \n\n',computationalTime);
        fprintf('_________________Computing The Error Distribution Ended__________________\n');
        fprintf('#########################################################################\n\n\n');
    end
end

end