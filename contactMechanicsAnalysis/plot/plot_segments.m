function [] = plot_segments(segments)
%   Plot segments on current graph
 
%% 
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
        % plot contours below the line
        plot(contourWall(:,1),contourWall(:,2),':','LineWidth',2,'Color','black');
        
        hold off;
    end
end

end