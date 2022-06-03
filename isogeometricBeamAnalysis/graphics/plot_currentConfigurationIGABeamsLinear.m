function index = plot_currentConfigurationIGABeamsLinear...
    (analysis, p, Xi, CP, isNURBS, dHat, rb, parameters, graph, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
% Plots the geometry together with the boundary conditions for a beam like 
% isogeometric structure together with the selected resultant for the
% current configuration
%
%      Input : 
%   analysis : Beam analysis type :
%                     'Bernoulli' : isogeometric Bernoulli beam analysis
%                    'Timoshenko' : isogeometric Timoshenko beam analysis
%        p,q : The polynomial degrees in u-,v-directions
%        U,V : The knot vectors in u-,v- directions
%         CP : The set of the Control points and weights
%    isNURBS : Flag on whether the geometrical basis is a NURBS or a 
%              B-Spline
%       dHat : The given displacement field
%         rb : The set of boundary conditions
% parameters : The technical parameters of the beam
%      graph : On the graphics
%     outMsg : Whether or not to output message on refinement progress
%              'outputEnabled' : enables output information
%
%  Output : 
%   index : The index of the current graph
%
% Function layout :
%
% 0. Read input
%
% 1. Create the first window (Deformed and/or undeformed geometry)
%
% 2. Create the second window (Resultant visualization)
%
% 3. Appendix
%
%% Function main body

if strcmp(outMsg,'outputEnabled')
    fprintf('____________________________________________________________________\n');
    fprintf('####################################################################\n');
    fprintf('Plotting the current configuration and the postprocessing resultants\n'); 
    fprintf('for the isogeometric ');
    if strcmp(analysis.type,'Bernoulli')
        fprintf('Bernoulli ');
    elseif strcmp(analysis.type,'Timoshenko')
        fprintf('Timoshenko ');
    end
    fprintf('beam analysis\n\n');
    if strcmp(graph.resultant,'displacement') && strcmp(graph.component,'x')
        fprintf('Chosen resultant to be displayed: Displacement\n');
        fprintf('Chosen component to be displayed: x-Component\n');
    elseif strcmp(graph.resultant,'displacement') && strcmp(graph.component,'y')
        fprintf('Chosen resultant to be displayed: Displacement\n');
        fprintf('Chosen component to be displayed: y-Component\n');
    elseif strcmp(graph.resultant,'displacement') && strcmp(graph.component,'2norm')
        fprintf('Chosen resultant to be displayed: Displacement\n');
        fprintf('Chosen component to be displayed: ||u||_2\n');
    elseif strcmp(graph.resultant,'crossSectionalRotation')
        fprintf('Chosen resultant to be displayed: Cross sectional rotation\n');
    elseif strcmp(graph.resultant,'force')
        fprintf('Chosen resultant to be displayed: Normal force\n');
    elseif strcmp(graph.resultant,'moment')
        fprintf('Chosen resultant to be displayed: Bending moment\n');
    elseif strcmp(graph.resultant,'shearForce')
        fprintf('Chosen resultant to be displayed: Shear force\n');
    end
    fprintf('____________________________________________________________________\n\n');

    % start measuring computational time
    tic;
end

%% 0. Read input

% Initialize figure
figure(graph.index);

% Assign the number of sampling points to be used for the postprocessing
gridPoints = 49;

% Number of evaluation points for the postprocessing resultants
numEval = 500;

% Initialize the array of the resultants
resultant = zeros(numEval+1,1);

%% 1. Create the first window (Deformed and/or undeformed geometry)

% Initialize subfigure 1
subplot(2,1,1);

% Compute the displaced Control Points
CPd = computeDisplacedControlPointsIGABeams2D(CP,dHat,analysis);

% Create arrays needed for the plotting of the surface, the supports as
% well as the force arrows
if strcmp(graph.postProcVisType,'undeformedShape')||strcmp(graph.postProcVisType,'deformedAndUndeformedShape')
    [Xp,Yp,Zp] = createBSplineCurveOnCartesianSpace(p,Xi,CP,isNURBS,gridPoints);
    if strcmp(analysis.type,'Bernoulli')
        [xs,ys,zs] = createSupportsForIGABernoulliBeam2D(CP,rb);
    elseif strcmp(analysis.type,'Timoshenko')
        [xs,ys,zs] = createSupportsForIGATimoshenkoBeam2D(CP,rb);
    end
end
if strcmp(graph.postProcVisType,'deformedShape')||strcmp(graph.postProcVisType,'deformedAndUndeformedShape')
    [XpDef,YpDef,ZpDef] = createBSplineCurveOnCartesianSpace(p,Xi,CPd,isNURBS,gridPoints);
    if ~strcmp(graph.postProcVisType,'deformedAndUndeformedShape')
        if strcmp(analysis.type,'Bernoulli')
            [xsDef,ysDef,zsDef] = createSupportsForIGABernoulliBeam2D(CPd,rb);
        elseif strcmp(analysis.type,'Timoshenko')
            [xsDef,ysDef,zsDef] = createSupportsForIGATimoshenkoBeam2D(CPd,rb);
        end
    end
end

% Plot the geometry, the knots on the geometry, the supports and the Control
% Polygon
if strcmp(graph.postProcVisType,'undeformedShape')||strcmp(graph.postProcVisType,'deformedAndUndeformedShape')
    % Plot geometry
    plot3(Xp,Yp,Zp,'Color','black');
    hold on;

    % Plot knots on the geometry
    plot_knotsForBSplineCurveOnCartesianSpace(p,Xi,CP,isNURBS);

    % Plot supports
    for k =1:length(xs(:,1))
        plot3(xs(k,:),ys(k,:),zs(k,:),'Linewidth',2,'Color','black');
    end

    % Plot the Control Points and the Control Polygon
    plot_ControlPolygonBSplineCurve(CP);
    
    if strcmp(graph.postProcVisType,'deformedAndUndeformedShape')
       hold on; 
    end
end
if strcmp(graph.postProcVisType,'deformedShape')||strcmp(graph.postProcVisType,'deformedAndUndeformedShape')
   % Plot geometry
    plot3(XpDef,YpDef,ZpDef,'Color','black');
    hold on;

    % Plot knots on the geometry
    plot_knotsForBSplineCurveOnCartesianSpace(p,Xi,CPd,isNURBS);

    % Plot supports
    if ~strcmp(graph.postProcVisType,'deformedAndUndeformedShape')
        for k =1:length(xsDef(:,1))
            plot3(xsDef(k,:),ysDef(k,:),zsDef(k,:),'Linewidth',2,'Color','black');
        end
    end

    % Plot the Control Points and the Control Polygon
    plot_ControlPolygonBSplineCurve(CPd);
end

% Choose plotting properties
% axis equal;
view(2);
camlight left; lighting phong;
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
zlabel('z','FontSize',14);
hold off;

% Choose a title for the subfigure
if strcmp(analysis.type,'Bernoulli')
    title('Deformed configuration for the isogeometric Bernoulli beam');
elseif strcmp(analysis.type,'Timoshenko')
    title('Deformed configuration for the isogeometric Timoshenko beam');
end

%% 2. Create the second window (Resultant visualization)

% Initialize subfigure 2
subplot(2,1,2);

% Compute a step
dxi = (Xi(length(Xi))-Xi(1))/numEval;

% Initialize the array of the coordinates on the parameter space
xi = zeros(numEval+1,1);

% Get the first coordinate on the NURBS parameter space
xi(1) = Xi(1);

% Loop over all sampling points
for j=1:numEval+1
    
    % Initializations
    dhatX = 0;
    dhatY = 0;
    crossSectionalRotation = 0;

    % Get a coordinate on the NURBS parameter space 
    if j ~= 1
        xi(j) = xi(j-1) + dxi; 
    end
    
    % Find the correct knot span
    knotSpan = findKnotSpan(xi(j),Xi,length(CP(:,1)));
    
    % Compute the NURBS basis functions and their derivatives
    if strcmp(graph.resultant,'displacement')
        noDeriv = 1;
    else
        noDeriv = 2;
    end
    dR = computeIGABasisFunctionsAndDerivativesForCurve...
        (knotSpan,p,xi(j),Xi,CP,isNURBS,noDeriv);

    if strcmp(graph.resultant,'displacement')
        % Compute iteratively the displacement field
        for b = 0:p
            if strcmp(analysis.type,'Bernoulli')
                dhatX = dR(b+1,1)*dHat(2*(knotSpan-p+b)-1,1) + dhatX;  
                dhatY = dR(b+1,1)*dHat(2*(knotSpan-p+b),1) + dhatY;
            elseif strcmp(analysis.type,'Timoshenko')
                dhatX = dR(b+1,1)*dHat(3*(knotSpan-p+b)-2,1) + dhatX;  
                dhatY = dR(b+1,1)*dHat(3*(knotSpan-p+b)-1,1) + dhatY;
            end
        end
        if strcmp(graph.component,'x')
            resultant(j) = dhatX;
        elseif strcmp(graph.component,'y')
            resultant(j) = dhatY;
        elseif strcmp(graph.component,'2norm')
            resultant(j) = sqrt(dhatX^2+dhatY^2);
        end
    elseif strcmp(graph.resultant,'crossSectionalRotation')
        if strcmp(analysis.type,'Bernoulli')
            error('The Bernoulli beam theory assumes no cross sectional deformation\n');
        elseif strcmp(analysis.type,'Timoshenko')
            for b = 0:p
                crossSectionalRotation = crossSectionalRotation + dR(b+1,1)*dHat(3*(knotSpan-p+b),1);
            end
            resultant(j) = crossSectionalRotation;
        end
    elseif strcmp(graph.resultant,'force')
        % Compute the base vector of the reference configuration
        [G,~] = computeBaseVectorAndNormalToNURBSCurve2D(knotSpan,p,xi,Xi,CP,dR);
        
        % Compute the normal strain at the given parametric location
        resultant(j) = computeNormalForceForIGABeams2D(knotSpan,p,dHat,parameters,analysis,G,dR);
    elseif strcmp(graph.resultant,'moment')
        % Compute the base vector, its derivative and the normal to the
        % curve and the acceleration vectors
        [G,dG] = computeBaseVectorNormalToNURBSCurveAndDeivativesForCurve2D(knotSpan,p,CP,dR);
        
        % Compute the bending moment at the parametric location
        resultant(j) = computeBendingMomentForIGABeams2D(knotSpan,p,dHat,parameters,analysis,G,dG,dR);
    elseif strcmp(graph.resultant,'shearForce')
        % Compute the base vector of the curve
        [G,~] = computeBaseVectorAndNormalToNURBSCurve2D(knotSpan,p,CPd,dR);
        
        % Compute the shear force at the parametric location
        resultant(j) = computeShearForceForIGABeams2D(knotSpan,p,dHat,parameters,analysis,G,dR);
    end
end

% Plot the resultant array over the parameter space
plot(xi,resultant);

% Choose graphic options
grid on;

% Adjust the title of the subfigure
if strcmp(graph.resultant,'displacement') && strcmp(graph.component,'x')
    title('Horizontal displacement u_x');
elseif strcmp(graph.resultant,'displacement') && strcmp(graph.component,'y')
    title('Horizontal displacement u_y');
elseif strcmp(graph.resultant,'displacement') && strcmp(graph.component,'2norm')
    title('Displacement magnitude ||u||_2');
elseif strcmp(graph.resultant,'crossSectionalRotation')
    title('Cross sectional rotation omega');
elseif strcmp(graph.resultant,'force')
    title('Normal force n');
elseif strcmp(graph.resultant,'moment')
    title('Bending moment m');
elseif strcmp(graph.resultant,'shearForce')
    title('Shear force q');
end

% Update the index counter
index = graph.index + 1;

%% 3. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;
    
    fprintf('Plotting the current configuration took %d seconds \n\n',computationalTime);
    fprintf('_________________Plotting Current Configuration Ended_______________\n');
    fprintf('####################################################################\n\n\n');
end

end
