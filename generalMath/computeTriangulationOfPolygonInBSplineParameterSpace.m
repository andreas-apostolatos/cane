function dt = computeTriangulationOfPolygonInBSplineParameterSpace ...
    (polygon, knotIntersections, parameterLines, tolPoints)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the triangulation resulting from the clipping of a polygonal are
% in the given 2D B-Spline parameter space.
%
%               Input :
%             polygon : The vertices of the polygon in a counterclockwise
%                       fashion polygon = zeros(noVertices,2)
%   knotIntersections : The knot intersections in the given B-Spline
%                       parameter space
%      parameterLines : The set of all the lines connecting the knots in
%                       the B-Spline parameter space
%           tolPoints : Tolerance for point coincidence while computing the
%                       intersections
%
%              Output :
%                  dt : Constrained Delaunay triangulation of the polygonal
%                       area,
%                                    .X : Vertices of the triangulation
%                        .Triangulation : The triangulation itself
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the connectivity lines in the B-Spline parameter space
% ->
%    1i. Loop over all the edges in the polygonal domain
%    ->
%        1i.1. Get the end points of the line in the polygon
%
%        1i.2. Get the end-points of the NURBS connectivity line
%
%        1i.3. Find possible intersection between the lines in the polygon and the B-Spline parameter space
%
%        1i.4. If intersection is found add it to the array
%    <-
% <-
%
% 2. If intersections have been found add them into the constraint vertices array
%
% 3. Add the polygon vertices into the constraint vertices array
%
% 4. Sort up the contraint vertices to a counterclockwise fashion
%
% 5. Get the knot intersections which are located inside the polygon
%
% 6. Add the constraint vertices into the point cloud to be processed by the Delaunay triangulation
%
% 7. Define the constraint edges by their end vertices [xIDLeft xIDRight] as they are sorted into the pointCloud array
%
% 8. Make the point cloud array unique but preserve the constraint edges in the same sequence at the beginning of the array
%
% 9. Triangulate the polygonal domain using a constraint Delaunay triangulation
%
%% Function main body

%% 0. Read input

% Initialize counter
counterIntersection = 1;

% Initialize the array of the intersections between the projected Finite
% Element and the knot lines
xiEtaIntersection = [];

% Initialize the bounding vertices of the polygonal area
boundingVertices = polygon;

%% 1. Check which knot intersections are in and one the projected polygon
[inPolygon,onPolygon] = inpolygon(knotIntersections(:,1),knotIntersections(:,2),polygon(:,1),polygon(:,2));
inPolygon(inPolygon == onPolygon) = false;
knotIntersectionsInPolygon = knotIntersections(inPolygon,:);
knotIntersectionsOnPolygon = knotIntersections(onPolygon,:);
if ~isempty(knotIntersectionsOnPolygon)
    boundingVertices = [boundingVertices
                        knotIntersectionsOnPolygon];
end

%% 2. Loop over all the connectivity lines in the B-Spline parameter space
for counterEdgesIGA = 1:length(parameterLines(:,1))
    %% 2i. Loop over all the edges in the polygonal domain
    for counterEdgesFEM = 1:length(polygon(:,1))
        %% 2i.1. Get the end points of the line in the polygon
        if counterEdgesFEM ~= length(polygon(:,1))
            xiEtaLeftFEM =  polygon(counterEdgesFEM,1:2);
            xiEtaRightFEM =  polygon(counterEdgesFEM + 1,1:2);
        else
            xiEtaLeftFEM =  polygon(length(polygon(:,1)),1:2);
            xiEtaRightFEM =  polygon(1,1:2);
        end
        
        %% 2i.2. Get the end-points of the NURBS connectivity line
        xiEtaLeftBSpline = parameterLines(counterEdgesIGA,1:2);
        xiEtaRightBSpline = parameterLines(counterEdgesIGA,3:4);
        
        %% 2i.3. Find possible intersection between the lines in the polygon and the B-Spline parameter space
        [xiEta,flag] = computeIntersectionBetweenStraightLines...
            (xiEtaLeftFEM,xiEtaRightFEM,xiEtaLeftBSpline,xiEtaRightBSpline);
        
        %% 2i.4. If intersection is found add it to the array
        if flag
            xiEtaIntersection(counterIntersection,:) = xiEta;
            counterIntersection = counterIntersection + 1;
        end
    end
end

%% 3. If intersections have been found add them into the constraint vertices array
if ~isempty(xiEtaIntersection)
    boundingVertices = [boundingVertices
                        xiEtaIntersection];
end

%% 4. Sort up the array of the bounding vertices in a counterclockwise fashion
boundingVertices = uniqueUpToTolerance(boundingVertices,tolPoints,'rows');
xi = boundingVertices(:,1);
eta = boundingVertices(:,2);
cxi = mean(xi);
ceta = mean(eta);
angles = atan2(eta - ceta, xi - cxi);
[~,I] = sort(angles);
boundingVertices = boundingVertices(I,:);

%% 5. Define the bounding vertices of the polygon
pointCloud = boundingVertices;
constraintEdges = zeros(length(pointCloud),2);
for counterConstraintEdges = 1:length(pointCloud)
    if counterConstraintEdges ~= length(pointCloud)
        constraintEdges(counterConstraintEdges,:) = ...
            [counterConstraintEdges counterConstraintEdges + 1];
    else
        constraintEdges(counterConstraintEdges,:) = ...
            [counterConstraintEdges 1];
    end
end

%% 6. Add the constraint vertices into the point cloud to be processed by the Delaunay triangulation
if ~isempty(knotIntersectionsInPolygon)
    pointCloud = [pointCloud
                  knotIntersectionsInPolygon];
end

%% 7. Triangulate the polygonal domain using a constraint Delaunay triangulation
dt = DelaunayTri(pointCloud,constraintEdges);

end
