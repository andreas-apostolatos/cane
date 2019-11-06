function index = plot_referenceConfigurationIGAThinStructureMultipatches...
    (BSplinePatches,connections,color,graph,outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function
%
% Plots the reference configuration for an isogeometric Kirchhoff-Love
% shell which includes the support triangles as well as the load arrows for
% a multipatch geometry.
%
%          Input : 
% BSplinePatches : Structure containing all the information regarding the 
%                  patches in the multipatch geometry
%    connections : Define the connection between the patches:
%                        .No : Number of connections
%                 .xiEtaCoup : [patchID1 patchID2 xi12 eta12 xi21 eta21
%                                 ...      ...    ...   ...  ...   ...]
%          graph : Structure containing information on the graphics
%                                  .index : The index of the current graph
%                     .isPrestressEnabled : Flag on whether the
%                                           visualization of the prestress
%                                           is on
%                          .compPrestress : Component of the prestress to
%                                           be visualized
%         outMsg : Enabled outputting information on the command window
%                  when chosen as 'outputEnabled'
%
%         Output : 
%          index : The index of the current graph
%         outMsg : Whether or not to output message on refinement progress
%                  'outputEnabled' : enables output information
%             
%                graphics
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the patches in the mesh
% ->
%    1i. Extract the patch information
%
%   1ii. Create the arrays containing the Cartesian coordinates of the B-Spline surfaces grid points, the supports and the arrow vectors
%
%  1iii. Plot the B-Spline surface
%
%   1iv. Plot the knots on the B-Spline surfaces
%
%    1v. Plot the control polygon of the B-Spline surfaces
%
%   1vi. Plot the supports on the geometries
%
%  1vii. Plot the load arrows on the geometries
% <-
%
% 2. Loop over all the connections
% ->
%    2i. Get the patch IDs
%
%   2ii. Get the coupling boundaries for both patches
%
%  2iii. Get the patch parameters
%
%   2iv. Compute a sequence of points over the interface curves
%
%    2v. Plot the coupling interface
% <-
%
% 3. Assign plotting properties
%
% 4. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('_____________________________________________________\n');
    fprintf('#####################################################\n');
    fprintf('Plotting the reference configuration for a multipatch\n');
    fprintf('isogeometric Kirchhoff-Love shell has been initiated\n');
    fprintf('_____________________________________________________\n\n');

    % start measuring computational time
    tic;
end

%% 0. Read input

% On the plotting of the element boundaries for the undeformed
% configuration
isDeformed = 0;

% Grid point number for the plotting of both the B-Spline surface and the
% knots
xiGrid = 29;
etaGrid = 29;

% Number of patches
noPatches = length(BSplinePatches);

% Initialize figure handle
figure(graph.index);

%% 1. Loop over all the patches in the mesh
for i = 1:noPatches
    %% 1i. Extract the patch information
    p = BSplinePatches{i}.p;
    q = BSplinePatches{i}.q;
    Xi = BSplinePatches{i}.Xi;
    Eta = BSplinePatches{i}.Eta;
    CP = BSplinePatches{i}.CP;
    isNURBS = BSplinePatches{i}.isNURBS;
    homDOFs = BSplinePatches{i}.homDOFs;
    inhomDOFs = BSplinePatches{i}.inhomDOFs;
    if isfield(graph,'isPrestressEnabled')
        if graph.isPrestressEnabled
            if isfield(BSplinePatches{i},'parameters')
                if isfield(BSplinePatches{i}.parameters,'prestress')
                    if ~isfield(graph,'compPrestress')
                        error('A component of the prestress for the visualization has to be defined');
                    end
                    prestress = BSplinePatches{i}.parameters.prestress;
                    compPrestress = graph.compPrestress;
                else
                    prestress = 'undefined';
                    compPrestress = 'undefined';
                end
            else
                prestress = 'undefined';
                compPrestress = 'undefined';
            end
        else
            prestress = 'undefined';
            compPrestress = 'undefined';
        end
    else
        prestress = 'undefined';
        compPrestress = 'undefined';
    end
    FGamma = BSplinePatches{i}.FGamma;

    %% 1ii. Create the arrays containing the Cartesian coordinates of the B-Spline surfaces grid points, the supports and the arrow vectors

    % Create the coordinates of the grid points on the surface
    [Xp,Yp,Zp,sigma0] = createBSplineSurfaceOnCartesianSpace...
        (p,q,Xi,Eta,CP,isNURBS,prestress,compPrestress,xiGrid,etaGrid);

    % Create the coordinates of the vertices of the support triangles for
    % the homogeneous Dirichlet boundary conditions
    [xsHom,ysHom,zsHom] = createSupports3D(CP,homDOFs);
    
    % Create the coordinates of the vertices of the support triangles for
    % the inhomogeneous Dirichlet boundary conditions
    [xsInhom,ysInhom,zsInhom] = createSupports3D(CP,inhomDOFs);
    
    % Create the start and end points of the arrows representing the loads
    if ischar(prestress)
        [xf,yf,zf] = createForceArrows3D(CP,FGamma);
    end

    %% 1iii. Plot the B-Spline surface
    if ischar(prestress)
        surf(Xp,Yp,Zp,'FaceColor',color,'EdgeColor','none');
    else
        surf(Xp,Yp,Zp,sigma0);
    end
    hold on;
    
    %% 1iv. Plot the knots on the B-Spline surfaces
    plot_knotsForBSplineSurfaceOnCartesianSpace(p,q,Xi,Eta,CP,isNURBS,isDeformed,xiGrid,etaGrid);
    
    %% 1v. Plot the control polygon of the B-Spline surfaces
%     if ischar(prestress)
%         plot_ControlPolygonBSplineSurface(CP);
%     end
    
    %% 1vi. Plot the supports on the geometries
    if ischar(prestress)
        for k = 1:length(xsHom(:,1))
            plot3(xsHom(k,:),ysHom(k,:),zsHom(k,:),'Linewidth',2,'Color','black');
        end
        for k = 1:length(xsInhom(:,1))
            plot3(xsInhom(k,:),ysInhom(k,:),zsInhom(k,:),'Linewidth',2,'Color','red');
        end
    end

    %% 1vii. Plot the load arrows on the geometries
    if ischar(prestress)
        for k = 1:length(xf(:,1))
            plot3(xf(k,:),yf(k,:),zf(k,:),'Color','blue','Linewidth',5);
            plot3(xf(k,1),yf(k,1),zf(k,1),'Marker','d','MarkerFaceColor','blue','MarkerSize',10);
        end
    end
end

%% 2. Loop over all the connections
if ~isempty(connections)
    for i = 1:connections.No
        %% 2i. Get the patch IDs

        % Patch 1 :
        % _________

        ID1 = connections.xiEtaCoup(i,1);

        % Patch 2 :
        % _________

        ID2 = connections.xiEtaCoup(i,2);

        %% 2ii. Get the coupling boundaries for both patches

        % Patch 1 :
        % _________

        xicoup1 = connections.xiEtaCoup(i,3:4);
        etacoup1 = connections.xiEtaCoup(i,5:6);

        % Patch 2 :
        % _________

        xicoup2 = connections.xiEtaCoup(i,7:8);
        etacoup2 = connections.xiEtaCoup(i,9:10);

        %% 2iii. Get the patch parameters

        % Patch 1 :
        % _________

        p1 = BSplinePatches{ID1}.p;
        q1 = BSplinePatches{ID1}.q;
        Xi1 = BSplinePatches{ID1}.Xi;
        Eta1 = BSplinePatches{ID1}.Eta;
        CP1 = BSplinePatches{ID1}.CP;
        isNURBS1 = BSplinePatches{ID1}.isNURBS;

        % Get the number of Control Points in each parametric direction
        nxi1 = length(CP1(:,1,1)); 
        neta1 = length(CP1(1,:,1));

        % Patch 2 :
        % _________

        p2 = BSplinePatches{ID2}.p;
        q2 = BSplinePatches{ID2}.q;
        Xi2 = BSplinePatches{ID2}.Xi;
        Eta2 = BSplinePatches{ID2}.Eta;
        CP2 = BSplinePatches{ID2}.CP;
        isNURBS2 = BSplinePatches{ID2}.isNURBS;

        % Get the number of Control Points in each parametric direction
        nxi2 = length(CP2(:,1,1)); 
        neta2 = length(CP2(1,:,1));

        %% 2iv. Compute a sequence of points over the interface curves

        % 1st patch :
        % ___________

        xi1 = xicoup1(1);
        eta1 = etacoup1(1);
        if etacoup1(1) == etacoup1(2)
            noPoints1 = xiGrid;
            dParameter1 = (xicoup1(2)-xicoup1(1))/noPoints1;
            etaSpan1 = findKnotSpan(eta1,Eta1,neta1);
        else
            noPoints1 = etaGrid;
            dParameter1 = (etacoup1(2)-etacoup1(1))/noPoints1;
            xiSpan1 = findKnotSpan(xi1,Xi1,nxi1);
        end

        % 2nd patch :
        % ___________

        xi2 = xicoup2(1);
        eta2 = etacoup2(1);
        if etacoup2(1) == etacoup2(2)
            nPoints2 = xiGrid;
            dParameter2 = (xicoup2(2)-xicoup2(1))/nPoints2;
            etaSpan2 = findKnotSpan(eta2,Eta2,neta2);
        else
            nPoints2 = etaGrid;
            dParameter2 = (etacoup2(2)-etacoup2(1))/nPoints2;
            xiSpan2 = findKnotSpan(xi2,Xi2,nxi2);
        end

        % Create the points on both curves

        % 1st patch :
        % ___________

        interfaceCurce1 = zeros(noPoints1,3);
        for j = 1:noPoints1
            % Find the span of the running parameter
            if etacoup1(1) == etacoup1(2)
                xiSpan1 = findKnotSpan(xi1,Xi1,nxi1);
            else
                etaSpan1 = findKnotSpan(eta1,Eta1,neta1);
            end

            % Compute the IGA basis functions at the surface parameters
            nDrv1 = 0;
            R1 = computeIGABasisFunctionsAndDerivativesForSurface...
                (xiSpan1,p1,xi1,Xi1,etaSpan1,q1,eta1,Eta1,CP1,isNURBS1,nDrv1);

            % Compute the point of the curve over the Cartesian space
            interfaceCurce1(j,:) = computeCartesianCoordinatesOfAPointOnBSplineSurface...
                (xiSpan1,p1,xi1,Xi1,etaSpan1,q1,eta1,Eta1,CP1,R1);

            % Update the running parameter
            if etacoup1(1) == etacoup1(2)
                xi1 = xi1 + dParameter1;
            else
                eta1 = eta1 + dParameter1;
            end
        end

        % 2nd patch :
        % ___________

        interfaceCurce2 = zeros(noPoints1,3);
        for j = 1:nPoints2
            % Find the span of the running parameter
            if etacoup2(1) == etacoup2(2)
                xiSpan2 = findKnotSpan(xi2,Xi2,nxi2);
            else
                etaSpan2 = findKnotSpan(eta2,Eta2,neta2);
            end

            % Compute the IGA basis functions at the surface parameters
            nDrv2 = 0;
            R2 = computeIGABasisFunctionsAndDerivativesForSurface...
                (xiSpan2,p2,xi2,Xi2,etaSpan2,q2,eta2,Eta2,CP2,isNURBS2,nDrv2);

            % Compute the point of the curve over the Cartesian space
            interfaceCurce2(j,:) = computeCartesianCoordinatesOfAPointOnBSplineSurface...
                (xiSpan2,p2,xi2,Xi2,etaSpan2,q2,eta2,Eta2,CP2,R2);

            % Update the running parameter
            if etacoup2(1) == etacoup2(2)
                xi2 = xi2 + dParameter2;
            else
                eta2 = eta2 + dParameter2;
            end
        end

        %% 2v. Plot the coupling interface

        % 1st patch :
        % ___________

        plot3(interfaceCurce1(:,1),interfaceCurce1(:,2),interfaceCurce1(:,3),...
            'LineWidth',2,'Color','magenta');

        % 2nd patch :
        % ___________

        plot3(interfaceCurce2(:,1),interfaceCurce2(:,2),interfaceCurce2(:,3),...
            'LineWidth',2,'Color','magenta');
    end
end

%% 3. Assign plotting properties
hold off;
axis equal;
if ~ischar(prestress)
    shading interp;
    colormap('jet');
    colorbar;
end
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
zlabel('z','FontSize',14);
index = graph.index + 1;

%% 4. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    fprintf('Plotting the reference configuration took %.2d seconds \n\n',computationalTime);
    fprintf('________Plotting Reference Configuration Ended_______\n');
    fprintf('#####################################################\n\n\n');
end

end
