function [] = plot_segments(segments)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Marko Leskovar
%                  Andreas Apostolatos
%
%% Function documentation
%
% Plot line segments on the current graph
%
%               Input : 
%     segments.points : Matrix with the coordinates of two wall determining
%                       points in form of [(x0,y0) ; (x1,y1)]. Length in
%                       direction 3 depends on the number of walls
%
%% Function main body

% check if segmets exist
if ~isempty(segments.points)
    
    % loop through the number of segments
    for i=1:size(segments.points,3)
        hold on;

        % plot a line form start to end coordiante x0 to x1
        startPoint = segments.points(:,1,i);
        endPoint = segments.points(:,2,i);
        plot(startPoint,endPoint,'-','Linewidth',2,'Color','black');
       
        % create contour
        contourWall = createWall(segments.points(:,:,i));
        % plot contours below the line to indicate orientation
        plot(contourWall(:,1),contourWall(:,2),':','LineWidth',2,'Color','black');
        
        hold off;
    end
end

end