function index = plot_referenceConfigurationIGABeams...
    (p,Xi,CP,isNURBS,rb,F,analysis,graph,outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Plots the geometry together with the boundary conditions and loads arrows
% for a beam like isogeometric structure of the Bernoulli or the Timoshenko
% type
%
%    Input : 
%        p : The polynomial degrees in u-,v-directions
%       Xi : The knot vectors in u-,v- directions
%       CP : The set of the Control points and weights
%  isNURBS : Flag on whether the basis is a NURBS or a B-Spline
%       rb : The set of boundary conditions
%        F : The force vector
% analysis : Beam analysis type :
%                   'Bernoulli' : isogeometric Bernoulli beam analysis
%                  'Timoshenko' : isogeometric Timoshenko beam analysis
%    graph : Structure on the graphics
%   outMsg : Whether or not to output message on refinement progress
%            'outputEnabled' : enables output information
%
%  Output : 
%   index : The index of the current graph
%
% Function layout :
%
% 0. Read input
%
% 1. Create the Cartesian coordinates of the points on the curve
%
% 2. Create the supports and the force arrows
%
% 3. Plot the reference geometry
%
% 4. Plot the knots on the geometry
%
% 5. Plot the supports on the geometry
%
% 6. Plot the load arrows on the geometry
%
% 7. Plot the Control Points and the Control Polygon
%
% 8. Define plotting properties
%
% 9. Update the index counter
%
% 10. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    if strcmp(analysis.type,'Bernoulli')
        fprintf('___________________________________________________________\n');
        fprintf('###########################################################\n');
    elseif strcmp(analysis.type,'Timoshenko')
        fprintf('____________________________________________________________\n');
        fprintf('############################################################\n');
    end
    fprintf('Plotting the reference configuration for the ');
    if strcmp(analysis.type,'Bernoulli')
        fprintf('Bernoulli ');
    elseif strcmp(analysis.type,'Timoshenko')
        fprintf('Timoshenko ');
    end
    fprintf('beam\n');
    if strcmp(analysis.type,'Bernoulli')
        fprintf('___________________________________________________________\n\n');
    elseif strcmp(analysis.type,'Timoshenko')
        fprintf('____________________________________________________________\n\n');
    end

    % start measuring computational time
    tic;
end

%% 0. Read input

% Assign a number of evaluation points
nEval = 49;

% Assign an index to the figure
figure(graph.index)

%% 1. Create the Cartesian coordinates of the points on the curve
[Xp,Yp,Zp] = createBSplineCurveOnCartesianSpace(p,Xi,CP,isNURBS,nEval);

%% 2. Create the supports and the force arrows
if strcmp(analysis.type,'Bernoulli')
    [xs,ys,zs] = createSupportsForIGABernoulliBeam2D(CP,rb);
    [xf,yf,zf] = createForceArrowsForIGABernoulliBeam2D(CP,F);
elseif strcmp(analysis.type,'Timoshenko')
    [xs,ys,zs] = createSupportsForIGATimoshenkoBeam2D(CP,rb);
    [xf,yf,zf] = createForceArrowsForIGATimoshenkoBeam2D(CP,F);
end

%% 3. Plot the reference geometry
plot3(Xp,Yp,Zp,'Color','black');
hold on;
 
%% 4. Plot the knots on the geometry
plot_knotsForBSplineCurveOnCartesianSpace(p,Xi,CP,isNURBS);
  
%% 5. Plot the supports on the geometry
for k =1:length(xs(:,1))
    plot3(xs(k,:),ys(k,:),zs(k,:),'Linewidth',2,'Color','black');
end
  
%% 6. Plot the load arrows on the geometry
for k =1:length(xf(:,1))
    plot3(xf(k,:),yf(k,:),zf(k,:),'Color','blue','Linewidth',5);
    plot3(xf(k,1),yf(k,1),zf(k,1),'Marker','d','MarkerFaceColor','blue','MarkerSize',10);
end

%% 7. Plot the Control Points and the Control Polygon
plot_ControlPolygonBSplineCurve(CP)

% Title of the figure
if strcmp(analysis.type,'Bernoulli')
    title('Reference configuration for an isogeometric Bernoulli beam');
elseif strcmp(analysis.type,'Timoshenko')
    title('Reference configuration for an isogeometric Timoshenko beam');
end

%% 8. Define plotting properties
axis equal;
view(2);
camlight left; lighting phong;
xlabel('X','FontSize',14);
ylabel('Y','FontSize',14);
zlabel('Z','FontSize',14);
hold off;

%% 9. Update the index counter
index = graph.index + 1;

%% 10. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;
    
    fprintf('Plotting the reference configuration took %d seconds \n\n',computationalTime);
    if strcmp(analysis.type,'Bernoulli')
        fprintf('___________Plotting Reference Configuration Ended__________\n');
        fprintf('###########################################################\n\n\n');
    elseif strcmp(analysis.type,'Timoshenko')
        fprintf('___________Plotting Reference Configuration Ended___________\n');
        fprintf('############################################################\n\n\n');
    end
end

end
